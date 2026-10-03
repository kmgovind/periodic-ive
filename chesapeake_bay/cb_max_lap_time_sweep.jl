# ==============================================================================
# CHESAPEAKE BAY MAXIMUM-LAP-TIME SWEEP
# ==============================================================================
# Standalone copy of the optimized dynamic Chesapeake Bay simulation.  This file
# deliberately does not include or modify cb_sim.jl.  It sweeps the maximum
# traversal duration in 24 h increments over the available CBOFS data window.
#
# Run from this directory:
#   julia --project=.. cb_max_lap_time_sweep.jl
# ==============================================================================

include("../src/cbdata.jl")
include("../src/ASVGeometry.jl")
include("../src/clarity.jl")
include("../src/IVESim.jl")

using .CBOFSData, .ASVGeometry, .Clarity, .IVESim
using Dates, Statistics, JuMP, Ipopt, JLD2, ProgressMeter, Plots, DelimitedFiles

# --- Chesapeake Bay configuration copied from cb_sim.jl ---
const U_NOMINAL = 1.75                 # m/s
const DT_STEP_SEC = 60.0               # s
const SIGMA_SENSING = 1500.0           # m
const ALPHA_CLARITY = 0.001
const P_IN_W_AVG = 750.0               # W
const WEIGHT_SIGMA = 2.0               # PSU
const SPEED_FLOOR = 0.0                # numerical guard; VehicleParams enforces 0.01 m/s
const SWEEP_INCREMENT_HOURS = 24
const OUTPUT_DIR = "max_lap_time_sweep"

"""Run the variable-speed dynamic simulation for one maximum lap duration."""
function run_dynamic_opt_sim_for_limit(frames, max_lap_time_sec)
    vehicle_params = VehicleParams(U_NOMINAL, 24 * 3600.0)

    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(
        transect_lon_shifted, transect_lat_shifted, npts
    )
    path_length_m = s_vec[end]
    n_segments = npts - 1
    ds = path_length_m / n_segments

    # Precompute local flat-Earth coordinates for the sensing calculation.
    r_earth = 6_371_000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    x_path = deg2rad.(lon_path) .* (r_earth * cos(mean_lat_rad))
    y_path = deg2rad.(lat_path) .* r_earth

    sim_times = collect(frames[1].dt:Minute(1):frames[end].dt)
    n_steps = length(sim_times)
    frame_times = [frame.dt for frame in frames]
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(frame_times .- t)))
    current_frame_idx = -1
    local lon_v, lat_v, salt_v

    clarity_history = zeros(Float64, npts, n_steps)
    weights_history = zeros(Float64, npts, n_steps)
    salinity_meas_hist = zeros(Float64, n_steps)
    lon_hist = zeros(Float64, n_steps)
    lat_hist = zeros(Float64, n_steps)
    lap_hist = zeros(Int, n_steps)
    speed_hist = zeros(Float64, n_steps)

    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2_000)
    sizehint!(lap_pos_buffer, 2_000)

    clarity_state = zeros(Float64, npts)
    current_weights = fill(1.0, npts)
    q_initial = fill(0.0, npts)
    u_initial = fill(U_NOMINAL, n_segments)
    s_robot_total = 0.0
    previous_lap = 1

    println("Solving first lap with T_lap <= ", max_lap_time_sec / 3600, " h")
    u_opt_segments = optimize_lap_speeds_full_spatiotemporal(
        current_weights, q_initial, n_segments, ds, P_IN_W_AVG, ALPHA_CLARITY,
        vehicle_params, SIGMA_SENSING;
        speed_floor=SPEED_FLOOR, max_lap_time_sec=max_lap_time_sec,
    )

    label = "T_max=$(round(Int, max_lap_time_sec / 3600)) h"
    @showprogress 1 "Simulating $label" for (t_idx, t) in enumerate(sim_times)
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m)

        if current_lap > previous_lap
            current_weights = calculate_target_weights(
                lap_sal_buffer, lap_pos_buffer, s_vec, npts; sigma_weight=WEIGHT_SIGMA,
            )
            u_initial = copy(u_opt_segments)
            u_opt_segments = optimize_lap_speeds_full_spatiotemporal_warm(
                current_weights, q_initial, u_initial, n_segments, ds, P_IN_W_AVG,
                ALPHA_CLARITY, vehicle_params, SIGMA_SENSING;
                speed_floor=SPEED_FLOOR, max_lap_time_sec=max_lap_time_sec,
            )
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            previous_lap = current_lap
        end

        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        f_idx = frame_idx_for_time(t)
        if f_idx != current_frame_idx
            lon_v, lat_v, salt_v = load_surface_salinity(frames[f_idx].file)
            current_frame_idx = f_idx
        end

        s_meas = sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)
        push!(lap_sal_buffer, s_meas)
        push!(lap_pos_buffer, s_mod)

        inner_dt = 0.25
        t_inner = 0.0
        s_t = s_robot_total
        u_current_display = 0.0
        while t_inner < DT_STEP_SEC
            dt_inc = min(inner_dt, DT_STEP_SEC - t_inner)
            s_t_mod = mod(s_t, path_length_m)
            seg_idx = min(floor(Int, s_t_mod / ds) + 1, n_segments)
            u_current = u_opt_segments[seg_idx]
            u_current_display = u_current

            x_t, y_t = get_2d_position(s_t_mod, s_vec, x_path, y_path)
            dist_sq = (x_path .- x_t).^2 .+ (y_path .- y_t).^2
            sensing = sensing_function.(dist_sq; S_0=0.005, sigma=SIGMA_SENSING)
            clarity_state .= update_clarity_rk4.(
                clarity_state, sensing, dt_inc; alpha=ALPHA_CLARITY,
            )
            s_t += u_current * dt_inc
            t_inner += dt_inc
        end

        s_robot_total = s_t
        clarity_history[:, t_idx] = clarity_state
        weights_history[:, t_idx] = current_weights
        salinity_meas_hist[t_idx] = s_meas
        lon_hist[t_idx] = lon_r
        lat_hist[t_idx] = lat_r
        lap_hist[t_idx] = current_lap
        speed_hist[t_idx] = u_current_display
    end

    return (; clarity_history, weights_history, salinity_meas_hist, lon_hist, lat_hist,
        lap_hist, speed_hist, sim_times, lon_path, lat_path, s_vec)
end

"""Reproduce the paper's speed/observed-weight plot for one sweep case."""
function plot_speed_weight_case(result, limit_hours, output_dir)
    n_steps = length(result.speed_hist)
    weight_at_robot = zeros(Float64, n_steps)
    for k in 1:n_steps
        d2 = (result.lon_path .- result.lon_hist[k]).^2 .+
             (result.lat_path .- result.lat_hist[k]).^2
        weight_at_robot[k] = result.weights_history[argmin(d2), k]
    end

    p_speed = plot(
        result.sim_times, result.speed_hist;
        ylabel="Speed (m/s)", title="T_lap <= $(limit_hours) h",
        grid=true, legend=false, ylims=(0, 1.1 * maximum(result.speed_hist)),
        linewidth=1.4,
    )
    p_weight = plot(
        result.sim_times, weight_at_robot;
        xlabel="Simulation time", ylabel="Weight at robot location",
        grid=true, legend=false, linewidth=1.4,
    )
    savefig(
        plot(p_speed, p_weight; layout=(2, 1), link=:x, size=(900, 650), dpi=300),
        joinpath(output_dir, "speed_weight_T$(limit_hours)h.pdf"),
    )
end

function main()
    mkpath(OUTPUT_DIR)
    frames = discover_files("./datafiles")
    isempty(frames) && error("No CBOFS NetCDF files found in ./datafiles")

    simulation_hours = floor(
        Int, Dates.value(frames[end].dt - frames[1].dt) / (1000 * 3600),
    )
    limits_hours = collect(SWEEP_INCREMENT_HOURS:SWEEP_INCREMENT_HOURS:simulation_hours)
    isempty(limits_hours) && error("Simulation window is shorter than 24 h")
    println("CBOFS window: $(simulation_hours) h; sweeping T_lap limits: $(limits_hours) h")

    speed_by_limit = Dict{Int, Vector{Float64}}()
    time_by_limit = Dict{Int, Vector{DateTime}}()
    summary_rows = Vector{Vector{Float64}}()

    for limit_hours in limits_hours
        result = run_dynamic_opt_sim_for_limit(frames, limit_hours * 3600.0)
        case_file = joinpath(OUTPUT_DIR, "T$(limit_hours)h_dynamic_opt.jld2")
        @save case_file result limit_hours

        plot_speed_weight_case(result, limit_hours, OUTPUT_DIR)
        speed_by_limit[limit_hours] = result.speed_hist
        time_by_limit[limit_hours] = result.sim_times

        total_clarity = sum(result.weights_history .* result.clarity_history)
        push!(summary_rows, [
            limit_hours,
            total_clarity,
            minimum(result.speed_hist),
            maximum(result.speed_hist),
            mean(result.speed_hist),
            maximum(result.lap_hist),
        ])
        println("T_lap <= $(limit_hours) h: total clarity=$(round(total_clarity; digits=2)), " *
                "speed range=$(round(minimum(result.speed_hist); digits=3))-" *
                "$(round(maximum(result.speed_hist); digits=3)) m/s, " *
                "laps=$(maximum(result.lap_hist))")
    end

    # One direct comparison of the speed trajectories across all lap-time limits.
    p_sweep = plot(
        xlabel="Simulation time", ylabel="Allocated speed (m/s)",
        title="Effect of maximum traversal time on allocated speed",
        grid=true, legend=:outerright, size=(1100, 550), dpi=300,
    )
    for limit_hours in limits_hours
        plot!(p_sweep, time_by_limit[limit_hours], speed_by_limit[limit_hours];
              label="T_lap <= $(limit_hours) h", linewidth=1.2)
    end
    savefig(p_sweep, joinpath(OUTPUT_DIR, "speed_profiles_by_max_lap_time.pdf"))

    header = ["max_lap_time_h", "total_clarity", "min_speed_mps", "max_speed_mps",
              "mean_speed_mps", "laps_completed"]
    writedlm(joinpath(OUTPUT_DIR, "summary.csv"), vcat(permutedims(header), reduce(vcat, permutedims.(summary_rows))), ',')
    println("Saved per-case data, paper-style plots, comparison plot, and summary to $(OUTPUT_DIR)/")
end

main()
