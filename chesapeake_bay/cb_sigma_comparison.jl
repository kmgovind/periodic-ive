# ==============================================================================
# CHESAPEAKE BAY TARGET-WEIGHT SIGMA COMPARISON
# ==============================================================================
# Standalone copy of the dynamic cb_sim workflow. It leaves cb_sim.jl unchanged.
# Each sigma uses a matched adaptive and constant-speed simulation. No maximum
# lap-time constraint is imposed; the 0.01 m/s VehicleParams bound is retained
# only to keep the travel-time parameterization finite.
#
# Run from this directory:
#   julia --project=.. cb_sigma_comparison.jl
# Regenerate the paper figures from the archived results without solving again:
#   PLOT_ONLY=true julia --project=.. cb_sigma_comparison.jl
# ==============================================================================

# Use the non-interactive GR backend so figure generation also works on a
# headless machine or continuous-integration runner.
ENV["GKSwstype"] = "100"

include("../src/cbdata.jl")
include("../src/ASVGeometry.jl")
include("../src/clarity.jl")
include("../src/IVESim.jl")

using .CBOFSData, .ASVGeometry, .Clarity, .IVESim
using Dates, Statistics, JuMP, Ipopt, JLD2, ProgressMeter, Plots, DelimitedFiles, Printf, LaTeXStrings

# Match the plotting style used throughout the paper.
const PLOT_FONT = "Computer Modern"
default(
    fontfamily=PLOT_FONT,
    linewidth=2,
    guidefontsize=14,
    tickfontsize=14,
    legendfontsize=10,
    titlefontsize=18,
)

# --- Chesapeake Bay setup, copied from cb_sim.jl ---
const U_NOMINAL = 1.75
const DT_STEP_SEC = 60.0
const SIGMA_SENSING = 1500.0
const ALPHA_CLARITY = 0.001
const P_IN_W_AVG = 750.0
const SPEED_FLOOR = 0.0       # VehicleParams supplies the 0.01 m/s numerical guard.
const SIGMA_VALUES = [2.0, 3.0, 4.0]
const OUTPUT_DIR = "sigma_comparison"
# Optional external export directory for a manuscript checkout. Leave unset for
# a standalone public-code clone; the canonical outputs always stay in
# OUTPUT_DIR. Example:
#   PAPER_FIGURE_DIR=/path/to/paper/figures julia --project=.. cb_sigma_comparison.jl
const PAPER_FIGURE_DIR = get(ENV, "PAPER_FIGURE_DIR", "")

function export_paper_figure(figure, filename)
    isempty(PAPER_FIGURE_DIR) && return
    mkpath(PAPER_FIGURE_DIR)
    savefig(figure, joinpath(PAPER_FIGURE_DIR, filename))
end

"""Shared geometry and time metadata for one dynamic simulation."""
function simulation_setup(frames)
    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(
        transect_lon_shifted, transect_lat_shifted, npts,
    )
    path_length_m = s_vec[end]
    n_segments = npts - 1
    ds = path_length_m / n_segments
    r_earth = 6_371_000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    x_path = deg2rad.(lon_path) .* (r_earth * cos(mean_lat_rad))
    y_path = deg2rad.(lat_path) .* r_earth
    sim_times = collect(frames[1].dt:Minute(1):frames[end].dt)
    frame_times = [frame.dt for frame in frames]
    return (; npts, lon_path, lat_path, s_vec, path_length_m, n_segments, ds,
        x_path, y_path, sim_times, frame_times)
end

"""Adaptive dynamic run for one target-weight sigma, with no lap-time cap."""
function run_dynamic_opt_sim(frames, sigma_weight)
    setup = simulation_setup(frames)
    vehicle_params = VehicleParams(U_NOMINAL, 24 * 3600.0)
    n_steps = length(setup.sim_times)
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(setup.frame_times .- t)))
    current_frame_idx = -1
    local lon_v, lat_v, salt_v

    clarity_state = zeros(Float64, setup.npts)
    clarity_history = zeros(Float64, setup.npts, n_steps)
    weights_history = zeros(Float64, setup.npts, n_steps)
    salinity_meas_hist = zeros(Float64, n_steps)
    lon_hist = zeros(Float64, n_steps)
    lat_hist = zeros(Float64, n_steps)
    lap_hist = zeros(Int, n_steps)
    speed_hist = zeros(Float64, n_steps)

    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2_000)
    sizehint!(lap_pos_buffer, 2_000)
    planned_profiles = Dict{Int, Vector{Float64}}()
    planned_weights = Dict{Int, Vector{Float64}}()

    current_weights = fill(1.0, setup.npts)
    q_initial = fill(0.0, setup.npts)
    u_initial = fill(U_NOMINAL, setup.n_segments)
    s_robot_total = 0.0
    previous_lap = 1

    println("Solving adaptive sigma=$(sigma_weight), lap 1 (no T_lap constraint)")
    u_opt_segments = optimize_lap_speeds_full_spatiotemporal(
        current_weights, q_initial, setup.n_segments, setup.ds, P_IN_W_AVG,
        ALPHA_CLARITY, vehicle_params, SIGMA_SENSING;
        speed_floor=SPEED_FLOOR, max_lap_time_sec=nothing,
    )
    planned_profiles[1] = copy(u_opt_segments)
    planned_weights[1] = copy(current_weights)

    @showprogress 1 "Adaptive sigma=$(sigma_weight)" for (t_idx, t) in enumerate(setup.sim_times)
        current_lap = floor(Int, s_robot_total / setup.path_length_m) + 1
        s_mod = mod(s_robot_total, setup.path_length_m)
        if current_lap > previous_lap
            current_weights = calculate_target_weights(
                lap_sal_buffer, lap_pos_buffer, setup.s_vec, setup.npts;
                sigma_weight=sigma_weight,
            )
            u_initial = copy(u_opt_segments)
            u_opt_segments = optimize_lap_speeds_full_spatiotemporal_warm(
                current_weights, q_initial, u_initial, setup.n_segments, setup.ds,
                P_IN_W_AVG, ALPHA_CLARITY, vehicle_params, SIGMA_SENSING;
                speed_floor=SPEED_FLOOR, max_lap_time_sec=nothing,
            )
            planned_profiles[current_lap] = copy(u_opt_segments)
            planned_weights[current_lap] = copy(current_weights)
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            previous_lap = current_lap
        end

        lon_r, lat_r = get_2d_position(s_mod, setup.s_vec, setup.lon_path, setup.lat_path)
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
            s_t_mod = mod(s_t, setup.path_length_m)
            seg_idx = min(floor(Int, s_t_mod / setup.ds) + 1, setup.n_segments)
            u_current = u_opt_segments[seg_idx]
            u_current_display = u_current
            x_t, y_t = get_2d_position(s_t_mod, setup.s_vec, setup.x_path, setup.y_path)
            dist_sq = (setup.x_path .- x_t).^2 .+ (setup.y_path .- y_t).^2
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

    return merge(setup, (; clarity_history, weights_history, salinity_meas_hist,
        lon_hist, lat_hist, lap_hist, speed_hist, planned_profiles, planned_weights))
end

"""Constant-speed dynamic baseline, evaluated with the same target-weight sigma."""
function run_dynamic_const_sim(frames, sigma_weight)
    setup = simulation_setup(frames)
    n_steps = length(setup.sim_times)
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(setup.frame_times .- t)))
    current_frame_idx = -1
    local lon_v, lat_v, salt_v

    clarity_state = zeros(Float64, setup.npts)
    clarity_history = zeros(Float64, setup.npts, n_steps)
    weights_history = zeros(Float64, setup.npts, n_steps)
    lap_hist = zeros(Int, n_steps)
    speed_hist = fill(U_NOMINAL, n_steps)
    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    s_robot_total = 0.0
    previous_lap = 1
    current_weights = fill(1.0, setup.npts)

    @showprogress 1 "Constant baseline sigma=$(sigma_weight)" for (t_idx, t) in enumerate(setup.sim_times)
        current_lap = floor(Int, s_robot_total / setup.path_length_m) + 1
        s_mod = mod(s_robot_total, setup.path_length_m)
        if current_lap > previous_lap
            current_weights = calculate_target_weights(
                lap_sal_buffer, lap_pos_buffer, setup.s_vec, setup.npts;
                sigma_weight=sigma_weight,
            )
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            previous_lap = current_lap
        end

        lon_r, lat_r = get_2d_position(s_mod, setup.s_vec, setup.lon_path, setup.lat_path)
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
        while t_inner < DT_STEP_SEC
            dt_inc = min(inner_dt, DT_STEP_SEC - t_inner)
            s_t_mod = mod(s_t, setup.path_length_m)
            x_t, y_t = get_2d_position(s_t_mod, setup.s_vec, setup.x_path, setup.y_path)
            dist_sq = (setup.x_path .- x_t).^2 .+ (setup.y_path .- y_t).^2
            sensing = sensing_function.(dist_sq; S_0=0.005, sigma=SIGMA_SENSING)
            clarity_state .= update_clarity_rk4.(
                clarity_state, sensing, dt_inc; alpha=ALPHA_CLARITY,
            )
            s_t += U_NOMINAL * dt_inc
            t_inner += dt_inc
        end

        s_robot_total = s_t
        clarity_history[:, t_idx] = clarity_state
        weights_history[:, t_idx] = current_weights
        lap_hist[t_idx] = current_lap
    end
    return merge(setup, (; clarity_history, weights_history, lap_hist, speed_hist))
end

function selected_completed_lap(result)
    # The last lap is typically incomplete at the end of the finite data window.
    return max(1, maximum(result.lap_hist) - 1)
end

"""Integrate the recorded speed history to cumulative path distance."""
function full_path_distance_km(result)
    distance_m = vcat(0.0, cumsum(result.speed_hist[1:end-1] .* DT_STEP_SEC))
    return distance_m ./ 1000
end

"""Plot all laps for each adaptive strategy on a common cumulative-distance axis."""
function plot_full_path_profiles(results, output_dir)
    ordered = sort(collect(results); by=first)
    max_distance_km = maximum(maximum(full_path_distance_km(result)) for (_, result) in ordered)
    panels = Any[]
    for (sigma_weight, result) in ordered
        p = plot(
            full_path_distance_km(result), result.speed_hist;
            xlabel="Cumulative distance along path (km)", ylabel="Speed (m/s)",
            title=latexstring("\\sigma_{\\zeta} = ", sigma_weight), grid=true, legend=false,
            xlims=(0, max_distance_km), ylims=(0, 4.6), color=:black,
            left_margin=10Plots.mm, bottom_margin=7Plots.mm,
            right_margin=4Plots.mm, top_margin=3Plots.mm,
        )
        push!(panels, p)
    end
    figure = plot(
        panels...; layout=(length(panels), 1), link=:x, size=(1100, 620), dpi=300,
        plot_title="Allocated speed over the full Chesapeake Bay trajectory",
    )
    savefig(figure, joinpath(output_dir, "speed_vs_full_path_by_sigma.pdf"))
    export_paper_figure(figure, "sigma_speed_vs_full_path.pdf")
end

"""Compare speed histories of all adaptive strategies and the constant baseline."""
function plot_speed_vs_time(adaptive_results, baseline_results, output_dir)
    p = plot(
        xlabel="Simulation time", ylabel="Speed (m/s)",
        title="Allocated speed over time for all strategies", grid=true,
        legend=:topright, size=(1100, 420), dpi=300,
        left_margin=10Plots.mm, bottom_margin=7Plots.mm,
        right_margin=5Plots.mm, top_margin=3Plots.mm,
    )
    for (sigma_weight, result) in sort(collect(adaptive_results); by=first)
        plot!(p, result.sim_times, result.speed_hist;
              label=latexstring("\\mathrm{IVE},\\ \\sigma_{\\zeta} = ", sigma_weight))
    end
    baseline_sigma = minimum(collect(keys(baseline_results)))
    baseline = baseline_results[baseline_sigma]
    plot!(p, baseline.sim_times, baseline.speed_hist;
          label="Constant-speed baseline", linestyle=:dash, color=:black)
    savefig(p, joinpath(output_dir, "speed_vs_time_by_strategy.pdf"))
    export_paper_figure(p, "sigma_speed_vs_time.pdf")
end

function write_paper_table(rows, output_dir)
    open(joinpath(output_dir, "sigma_comparison_table.tex"), "w") do io
        println(io, "% Replace the table body in simulation-study.tex after reviewing values.")
        println(io, "\\begin{tabular}{c|r|r|r|r|r}")
        println(io, "\\hline")
        println(io, "\$\\sigma_\\zeta\$ & Clarity (IVE) & Clarity (constant) & Gain & Min. speed & Lap time " * repeat(string(Char(92)), 2))
        println(io, "\\hline")
        for row in rows
            line = @sprintf("%.1f & %.2f & %.2f & %.1f\\%% & %.3f & %.2f h ",
                            row[1], row[2], row[3], row[4], row[5], row[6])
            println(io, line * repeat(string(Char(92)), 2))
        end
        println(io, "\\hline")
        println(io, "\\end{tabular}")
    end
end

function main()
    mkpath(OUTPUT_DIR)
    frames = discover_files("./datafiles")
    isempty(frames) && error("No CBOFS NetCDF files found in ./datafiles")

    adaptive_results = Dict{Float64, Any}()
    baseline_results = Dict{Float64, Any}()
    summary_rows = Vector{Vector{Float64}}()
    for sigma_weight in SIGMA_VALUES
        adaptive = run_dynamic_opt_sim(frames, sigma_weight)
        baseline = run_dynamic_const_sim(frames, sigma_weight)
        adaptive_results[sigma_weight] = adaptive
        baseline_results[sigma_weight] = baseline
        adaptive_file = joinpath(OUTPUT_DIR, "sigma_$(sigma_weight)_adaptive.jld2")
        baseline_file = joinpath(OUTPUT_DIR, "sigma_$(sigma_weight)_constant.jld2")
        @save adaptive_file adaptive sigma_weight
        @save baseline_file baseline sigma_weight

        lap = selected_completed_lap(adaptive)
        profile = adaptive.planned_profiles[lap]
        clarity_adaptive = sum(adaptive.weights_history .* adaptive.clarity_history)
        clarity_constant = sum(baseline.weights_history .* baseline.clarity_history)
        gain_percent = 100 * (clarity_adaptive - clarity_constant) / abs(clarity_constant)
        planned_lap_time_h = sum(adaptive.ds ./ profile) / 3600
        # Report full-run, rather than selected-lap, low-speed behavior.
        floor_active_fraction = mean(adaptive.speed_hist .<= 0.0101)
        push!(summary_rows, [sigma_weight, clarity_adaptive, clarity_constant, gain_percent,
                             minimum(adaptive.speed_hist), planned_lap_time_h, maximum(adaptive.speed_hist),
                             floor_active_fraction, lap])
        println(@sprintf("sigma=%.1f | gain=%.2f%% | full-run min u=%.3f m/s | T_lap=%.2f h | floor active=%.1f%%",
                         sigma_weight, gain_percent, minimum(adaptive.speed_hist), planned_lap_time_h,
                         100 * floor_active_fraction))
    end

    plot_full_path_profiles(adaptive_results, OUTPUT_DIR)
    plot_speed_vs_time(adaptive_results, baseline_results, OUTPUT_DIR)
    header = ["sigma_weight", "clarity_adaptive", "clarity_constant", "gain_percent",
              "min_speed_mps", "selected_lap_time_h", "max_speed_mps",
              "floor_active_fraction", "selected_completed_lap"]
    writedlm(joinpath(OUTPUT_DIR, "sigma_comparison_summary.csv"),
             vcat(permutedims(header), reduce(vcat, permutedims.(summary_rows))), ',')
    write_paper_table(summary_rows, OUTPUT_DIR)
    println("Saved results, position-profile plot, CSV summary, and LaTex table to $(OUTPUT_DIR)/")
end

"""Regenerate figures from completed simulations without rerunning the sweep."""
function plot_existing_results()
    adaptive_results = Dict{Float64, Any}()
    baseline_results = Dict{Float64, Any}()
    summary_rows = Vector{Vector{Float64}}()
    for sigma_weight in SIGMA_VALUES
        adaptive_file = joinpath(OUTPUT_DIR, "sigma_$(sigma_weight)_adaptive.jld2")
        baseline_file = joinpath(OUTPUT_DIR, "sigma_$(sigma_weight)_constant.jld2")
        isfile(adaptive_file) || error("Missing completed result: $(adaptive_file)")
        isfile(baseline_file) || error("Missing completed result: $(baseline_file)")
        adaptive_results[sigma_weight] = load(adaptive_file, "adaptive")
        baseline_results[sigma_weight] = load(baseline_file, "baseline")
        adaptive = adaptive_results[sigma_weight]
        baseline = baseline_results[sigma_weight]
        lap = selected_completed_lap(adaptive)
        profile = adaptive.planned_profiles[lap]
        clarity_adaptive = sum(adaptive.weights_history .* adaptive.clarity_history)
        clarity_constant = sum(baseline.weights_history .* baseline.clarity_history)
        gain_percent = 100 * (clarity_adaptive - clarity_constant) / abs(clarity_constant)
        push!(summary_rows, [sigma_weight, clarity_adaptive, clarity_constant, gain_percent,
                             minimum(adaptive.speed_hist), sum(adaptive.ds ./ profile) / 3600,
                             maximum(adaptive.speed_hist), mean(adaptive.speed_hist .<= 0.0101), lap])
    end
    plot_full_path_profiles(adaptive_results, OUTPUT_DIR)
    plot_speed_vs_time(adaptive_results, baseline_results, OUTPUT_DIR)
    header = ["sigma_weight", "clarity_adaptive", "clarity_constant", "gain_percent",
              "min_speed_mps", "selected_lap_time_h", "max_speed_mps",
              "floor_active_fraction", "selected_completed_lap"]
    writedlm(joinpath(OUTPUT_DIR, "sigma_comparison_summary.csv"),
             vcat(permutedims(header), reduce(vcat, permutedims.(summary_rows))), ',')
    println("Regenerated figures from completed results in $(OUTPUT_DIR)/")
end

if get(ENV, "PLOT_ONLY", "false") == "true"
    plot_existing_results()
else
    main()
end
