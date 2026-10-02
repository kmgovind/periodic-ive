# Importance-weight decay sensitivity for the Chesapeake Bay adaptive controller.
#
# Larger sigma_weight means a slower Gaussian decline away from the rated
# salinity and therefore more spatially even importance weights.  The study
# keeps every physical/model parameter fixed and asks whether the imposed
# 0.5 m/s floor remains active in each re-optimization.
#
# Run from this directory:
#   julia --project=../ importance_weight_decay_sensitivity.jl
# Optional: SENSITIVITY_HOURS=72 SIGMA_WEIGHTS=1.5,2.0,3.0,4.5 julia --project=../ importance_weight_decay_sensitivity.jl

ENV["GKSwstype"] = "100"  # non-interactive plots

include("../src/cbdata.jl")
include("../src/ASVGeometry.jl")
include("../src/clarity.jl")
include("../src/IVESim.jl")

using .CBOFSData, .ASVGeometry, .Clarity, .IVESim
using Dates, Statistics, DataFrames, JLD2, Plots, Printf

const SPEED_FLOOR = 0.5
const FLOOR_TOL = 1e-3
const U_NOMINAL = 1.75
const ALPHA_CLARITY = 0.001
const SIGMA_SENSING = 1500.0
const NPTS = 100                     # screening resolution; use 200 for a final confirmation
const DT_STEP = 60.0

parse_float_list(text) = parse.(Float64, split(text, ','))
sigma_weights = parse_float_list(get(ENV, "SIGMA_WEIGHTS", "1.5,2.0,3.0,4.5,6.0"))
max_hours = parse(Float64, get(ENV, "SENSITIVITY_HOURS", "48"))

"Summarize whether an optimized segment profile is constrained by the speed floor."
function floor_metrics(u, ds; floor=SPEED_FLOOR, tol=FLOOR_TOL)
    dt = ds ./ u
    active = u .<= floor + tol
    return (
        min_speed_mps = minimum(u),
        floor_active_segments = count(active),
        floor_active_fraction = mean(active),
        floor_active_time_fraction = sum(dt[active]) / sum(dt),
        lap_time_hr = sum(dt) / 3600,
    )
end

"Run the dynamic simulation for one salinity-weight width and log every optimization." 
function run_case(frames, sigma_weight; max_hours=max_hours)
    vehicle = VehicleParams(U_NOMINAL, 24 * 3600.0)
    P_in = vehicle.kh + vehicle.km * U_NOMINAL^3
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, NPTS)
    path_length = s_vec[end]
    N = NPTS - 1
    ds = path_length / N

    start_time = frames[1].dt
    end_time = min(frames[end].dt, start_time + Millisecond(round(Int, max_hours * 3600_000)))
    sim_times = collect(start_time:Minute(1):end_time)
    frame_times = [f.dt for f in frames]
    frame_index(t) = argmin(abs.(Dates.value.(frame_times .- t)))

    # Cartesian path coordinates used by the sensing model.
    earth_radius = 6_371_000.0
    scale_x = earth_radius * cos(deg2rad(mean(lat_path)))
    X = deg2rad.(lon_path) .* scale_x
    Y = deg2rad.(lat_path) .* earth_radius

    q = zeros(NPTS)
    weights = ones(NPTS)
    u = optimize_lap_speeds_full_spatiotemporal(
        weights, q, N, ds, P_in, ALPHA_CLARITY, vehicle, SIGMA_SENSING;
        speed_floor=SPEED_FLOOR,
    )

    rows = DataFrame(sigma_weight=Float64[], lap=Int[], kind=String[],
        weight_mean=Float64[], weight_std=Float64[], weight_cv=Float64[],
        min_speed_mps=Float64[], floor_active_segments=Int[],
        floor_active_fraction=Float64[], floor_active_time_fraction=Float64[], lap_time_hr=Float64[])
    profiles = DataFrame(sigma_weight=Float64[], lap=Int[], kind=String[], segment=Int[],
        distance_m=Float64[], speed_mps=Float64[], importance_weight=Float64[])
    function record!(lap, kind, profile)
        m = floor_metrics(profile, ds)
        push!(rows, (sigma_weight, lap, kind, mean(weights), std(weights), std(weights)/mean(weights),
            m.min_speed_mps, m.floor_active_segments, m.floor_active_fraction,
            m.floor_active_time_fraction, m.lap_time_hr))
        for i in eachindex(profile)
            push!(profiles, (sigma_weight, lap, kind, i, (i - 0.5) * ds, profile[i], weights[i]))
        end
    end
    record!(1, "uniform", u)

    position_total = 0.0
    previous_lap = 1
    salinity_buffer = Float64[]
    position_buffer = Float64[]
    cached_frame = -1
    lon_v = lat_v = salt_v = nothing

    for t in sim_times
        lap = floor(Int, position_total / path_length) + 1
        position_mod = mod(position_total, path_length)
        if lap > previous_lap
            weights = calculate_target_weights(salinity_buffer, position_buffer, s_vec, NPTS;
                sigma_weight=sigma_weight)
            u = optimize_lap_speeds_full_spatiotemporal_warm(
                weights, q, u, N, ds, P_in, ALPHA_CLARITY, vehicle, SIGMA_SENSING;
                speed_floor=SPEED_FLOOR,
            )
            record!(lap, "adaptive", u)
            empty!(salinity_buffer); empty!(position_buffer)
            previous_lap = lap
        end

        idx = frame_index(t)
        if idx != cached_frame
            lon_v, lat_v, salt_v = load_surface_salinity(frames[idx].file)
            cached_frame = idx
        end
        lon_r, lat_r = get_2d_position(position_mod, s_vec, lon_path, lat_path)
        push!(salinity_buffer, sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v))
        push!(position_buffer, position_mod)

        # Match the full simulation's 0.25-s sensing/propagation update.
        elapsed = 0.0
        position_local = position_total
        while elapsed < DT_STEP
            h = min(0.25, DT_STEP - elapsed)
            s_local = mod(position_local, path_length)
            segment = min(floor(Int, s_local / ds) + 1, N)
            Xr, Yr = get_2d_position(s_local, s_vec, X, Y)
            sensing = sensing_function.((X .- Xr).^2 .+ (Y .- Yr).^2; S_0=0.005, sigma=SIGMA_SENSING)
            q .= update_clarity_rk4.(q, sensing, h; alpha=ALPHA_CLARITY)
            position_local += u[segment] * h
            elapsed += h
        end
        position_total = position_local
    end
    return (; metrics=rows, profiles)
end

frames = discover_files("./datafiles")
@info "Importance-weight sensitivity" max_hours sigma_weights speed_floor=SPEED_FLOOR npts=NPTS
all_rows = DataFrame()
all_profiles = DataFrame()
for sigma_weight in sigma_weights
    @info "Running case" sigma_weight
    result = run_case(frames, sigma_weight)
    append!(all_rows, result.metrics)
    append!(all_profiles, result.profiles)
end

# Adaptive rows are the direct answer to the PI's question. Aggregate across
# re-optimizations rather than across minute-by-minute samples of a profile.
adaptive = filter(:kind => ==("adaptive"), all_rows)
# A short smoke test may not finish its initial lap. Preserve that result instead
# of failing during aggregation, but label it so it is not mistaken for evidence
# about adaptive importance weights.
analysis_rows = isempty(adaptive) ? all_rows : adaptive
analysis_kind = isempty(adaptive) ? "initial uniform optimization only" : "adaptive re-optimizations"
isempty(adaptive) && @warn "No lap completed in this run; increase SENSITIVITY_HOURS before drawing conclusions about sigma_weight."
summary = combine(groupby(analysis_rows, :sigma_weight),
    nrow => :adaptive_reoptimizations,
    :floor_active_fraction => mean => :mean_floor_active_fraction,
    :floor_active_fraction => maximum => :max_floor_active_fraction,
    :floor_active_time_fraction => mean => :mean_floor_active_time_fraction,
    :floor_active_segments => mean => :mean_floor_active_segments,
    :min_speed_mps => minimum => :minimum_speed_seen_mps,
    :weight_cv => mean => :mean_weight_cv,
)
sort!(summary, :sigma_weight)
println("\n=== Adaptive speed-floor sensitivity summary ===")
show(summary, allrows=true, allcols=true); println()

stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
@save "importance_weight_decay_sensitivity_$stamp.jld2" all_rows all_profiles summary sigma_weights max_hours

p_floor = plot(analysis_rows.sigma_weight, analysis_rows.floor_active_fraction, group=analysis_rows.lap,
    seriestype=:scatter, xlabel="Weight width σ_weight (larger = slower decay)",
    ylabel="Fraction of segments at 0.5 m/s", label="lap")
p_time = plot(summary.sigma_weight, summary.mean_floor_active_time_fraction, marker=:circle,
    xlabel="Weight width σ_weight (larger = slower decay)", ylabel="Mean time fraction at speed floor",
    label="adaptive mean")
p_cv = plot(summary.sigma_weight, summary.mean_weight_cv, marker=:circle,
    xlabel="Weight width σ_weight (larger = slower decay)", ylabel="Mean coefficient of variation of weights",
    label="weight heterogeneity")
figure = plot(p_floor, p_time, p_cv, layout=(3,1), size=(850,950))
savefig(figure, "importance_weight_decay_sensitivity_$stamp.png")

# One directly interpretable speed-allocation figure: the final adaptive profile
# for each decay setting, with the associated importance weights beneath it.
p_speed = plot(xlabel="Distance along transect (m)", ylabel="Allocated speed (m/s)",
    title="Final adaptive speed allocation", legend=:bottomright)
p_weight = plot(xlabel="Distance along transect (m)", ylabel="Importance weight",
    title="Associated target-weight profile", legend=:bottomright)
for σw in sigma_weights
    available = filter(row -> row.sigma_weight == σw && row.kind == "adaptive", all_profiles)
    isempty(available) && continue
    last_lap = maximum(available.lap)
    profile = filter(row -> row.lap == last_lap, available)
    plot!(p_speed, profile.distance_m, profile.speed_mps, label="σ=$(σw)")
    plot!(p_weight, profile.distance_m, profile.importance_weight, label="σ=$(σw)")
end
savefig(plot(p_speed, p_weight, layout=(2,1), size=(950,700)),
    "importance_weight_decay_speed_profiles_$stamp.png")

println("\nSaved raw per-lap metrics, summary, and figure with timestamp $stamp.")
println("Aggregation used: $analysis_kind.")
println("Interpretation: for adaptive re-optimizations, a lower active fraction/time fraction as σ_weight rises supports the PI's hypothesis.")
