if !isdefined(@__MODULE__, :CBOFSData)
    include("../src/cbdata.jl")
    include("../src/ASVGeometry.jl")
    include("../src/clarity.jl")
    include("../src/IVESim.jl")
    include("../src/ntfy.jl")
    using .CBOFSData, .ASVGeometry, .Clarity, .IVESim, .Notify
    
    using Dates, Statistics, JuMP, Ipopt, JLD2, ProgressMeter, Base.Threads, Printf, Plots
end

u_nominal = 1.75 # m/s
dt_step_sec = 60.0 # seconds
sigma_val = 300.0 # sensing footprint in meters
alpha_val = 0.001 # clarity decay rate

function run_cb_sim()
    # ==============================================================================
    # 1. INITIALIZATION & SETUP
    # ==============================================================================
    vehicle_params = VehicleParams(u_nominal, 24*3600.0)

    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]

    # Precompute Local Flat-Earth Cartesian Coordinates for the Path
    R_earth = 6371000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    cos_mean_lat = cos(mean_lat_rad)

    X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
    Y_path = deg2rad.(lat_path) .* R_earth

    power_nominal = vehicle_params.kh + vehicle_params.km * (u_nominal^3)
    P_in_W = power_nominal          # Continuous Solar Power Generation [Watts]

    # Path Discretization for JuMP
    N_segments = npts - 1
    ds = path_length_m / N_segments

    # Generate times 
    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)

    # Initial clarity state (The physical transient starts at 0)
    clarity_state = zeros(Float64, npts)

    # History Arrays
    clarity_history    = zeros(Float64, npts, N_steps)
    weights_history    = zeros(Float64, npts, N_steps)
    salinity_meas_hist = zeros(Float64, N_steps)
    lon_hist           = zeros(Float64, N_steps)
    lat_hist           = zeros(Float64, N_steps)
    lap_hist           = zeros(Int, N_steps)
    speed_hist         = zeros(Float64, N_steps)

    # Buffers for the current lap's measurements
    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2000) 
    sizehint!(lap_pos_buffer, 2000)

    s_robot_total = 0.0 
    previous_lap = 1 
                      
    # For static environment: always use the first frame
    frame_idx_for_time(t) = 1

    # ==============================================================================
    # 2. LAP 1 OFFLINE OPTIMIZATION (Paradigm 1)
    # ==============================================================================
    current_weights = fill(1.0, npts)
    q_initial = fill(0.0, npts)  # Initial clarity guess for optimization
    u_initial = fill(u_nominal, N_segments) # Initial speed guess for warm start
    opt_time_start = time()
    
    # Notice: We do NOT pass clarity_state. The solver just finds the periodic orbit!
    u_opt_segments = optimize_lap_speeds_full_spatiotemporal_warm(
        current_weights, q_initial, u_initial, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val
    )
    
    opt_time_elapsed = time() - opt_time_start
    println("Initial optimization completed in $(round(opt_time_elapsed, digits=2)) seconds")

    # ==============================================================================
    # 3. MAIN SIMULATION LOOP
    # ==============================================================================
    println("Starting Optimized Variable-Speed Simulation Loop...")

    @showprogress 1 "Running Sim..." for (t_idx, t) in enumerate(sim_times)
        
        # --- A. Determine Lap & Position ---
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m) 
        
        # --- B. LAP TRANSITION LOGIC (Adaptive Paradigm 1) ---
        if current_lap > previous_lap
            # 1. Update weights based on collected lap data
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)

            # 2. WARM START: Use the previous lap's optimal speeds as the new guess
            u_initial = copy(u_opt_segments)

            # 3. Re-solve for the NEW infinite-horizon limit cycle
            u_opt_segments = optimize_lap_speeds_full_spatiotemporal_warm(
                current_weights, q_initial, u_initial, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val
            )
            
            # 4. Clear buffers for the next lap
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            
            previous_lap = current_lap
        end
        
        # --- C. Environmental Lookup & Buffering ---
        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        # Use static salinity arrays loaded once at startup
        lon_v, lat_v, salt_v = lon_v_static, lat_v_static, salt_v_static
        s_meas = sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)
        
        push!(lap_sal_buffer, s_meas)
        push!(lap_pos_buffer, s_mod)

        # =====================================================================
        # --- D. HIGH-RESOLUTION PHYSICS SWEEP ---
        # =====================================================================
        # Even though the optimizer assumes periodic steady-state, our physics 
        # engine rigidly simulates the continuous transient behavior!
        inner_dt = 0.25 
        t_inner = 0.0
        path_length_true = s_vec[end]
        s_t = s_robot_total
        u_current_display = 0.0
        
        while t_inner < dt_step_sec
            dt_inc = min(inner_dt, dt_step_sec - t_inner)
             
            s_t_mod = mod(s_t, path_length_true)
            
            seg_idx = min(floor(Int, s_t_mod / ds) + 1, N_segments)
            u_current = u_opt_segments[seg_idx]
            u_current_display = u_current
            
            X_t, Y_t = get_2d_position(s_t_mod, s_vec, X_path, Y_path)
            
            dx = X_path .- X_t
            dy = Y_path .- Y_t
            dist_sq = dx.^2 .+ dy.^2
            
            Sj_current = sensing_function.(dist_sq; S_0 = 1.0, sigma = sigma_val)
            clarity_state .= update_clarity_rk4.(clarity_state, Sj_current, dt_inc; alpha=alpha_val)
            
            s_t += u_current * dt_inc
            t_inner += dt_inc
        end
        
        s_robot_total = s_t
        
        # --- E. Log Everything ---
        clarity_history[:, t_idx]  = clarity_state
        weights_history[:, t_idx]  = current_weights
        salinity_meas_hist[t_idx]  = s_meas
        lon_hist[t_idx]            = lon_r
        lat_hist[t_idx]            = lat_r
        lap_hist[t_idx]            = current_lap
        speed_hist[t_idx]          = u_current_display
    end

    return clarity_history, speed_hist, lap_hist, lon_hist, lat_hist, salinity_meas_hist, weights_history
end
################### Transect Definition #######################################################################################
### Load salinity data
frames = discover_files("./datafiles")

# Use a static environment: load the first frame once and reuse for all timesteps
first_frame = frames[25]
lon_v_static, lat_v_static, salt_v_static = load_surface_salinity(first_frame.file)

fr = frames[1]
lon_valid, lat_valid, salt_valid = load_surface_salinity(fr.file)
roi = (lat_valid .>= lat_min) .& (lat_valid .<= lat_max)

# Run variable speed simulation
try
    println("Running simulation...")
    clarity_history, speed_hist, lap_hist, lon_hist, lat_hist, salinity_meas_hist, weights_history = run_cb_sim()

    # Save results to JLD2 file for later analysis
    filename = Dates.format(now(), "yyyy-mm-dd_HHMMSS") * "_static_sim_results.jld2"
    @save filename clarity_history speed_hist lap_hist lon_hist lat_hist salinity_meas_hist weights_history
    println("Simulation results saved to ", filename)


    npts = size(clarity_history, 1)
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    frame_times = [fr.dt for fr in frames]
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(frame_times .- t)))

    ################################### Plotting & Analysis ########################################
    total_clarity = sum(weights_history .* clarity_history)
    println("Total Accumulated Clarity (All Time): ", round(total_clarity, digits=2))

    total_laps = lap_hist[end]
    lap_clarity = zeros(Float64, total_laps)

    for lap in 1:total_laps
        idx_for_lap = findall(x -> x == lap, lap_hist)
        lap_clarity[lap] = sum(clarity_history[:, idx_for_lap])
        println("Lap $lap Total Clarity: ", round(lap_clarity[lap], digits=2))
    end

    target_salinity = 30.0
    salinity_deviation = salinity_meas_hist .- target_salinity

    println("Mean Salinity Deviation: ", round(mean(salinity_deviation), digits=2), " PSU")
    println("Standard Deviation of Salinity Deviation: ", round(std(salinity_deviation), digits=2), " PSU")

    histogram(
        salinity_deviation,
        bins = 30,
        normalize = :probability,
        title = "Deviation from Target Salinity (30 PSU)\nStatic Environment",
        xlabel = "Deviation (PSU)",
        ylabel = "Frequency",
        label = "Measurements",
        color = :steelblue,
        linecolor = :white,
        grid = true,
        size = (600, 400)
    )
    vline!([0.0], lw=3, color=:red, label="Target (30 PSU)")
    savefig("static_env_speed_salinity_deviation_histogram.png")


    nsteps = length(sim_times)
    @assert nsteps > 0 "No simulation history found. Run the simulation cell first."

    frame_stride = max(1, fld(nsteps, 300))
    frame_ids = 1:frame_stride:nsteps

    lat_min, lat_max = 36.75, 37.5
    lon_pad, lat_pad = 0.02, 0.02
    xlims_map = (minimum(lon_path) - lon_pad, maximum(lon_path) + lon_pad)
    ylims_map = (lat_min - lat_pad, lat_max + lat_pad)

    smin = minimum(salinity_meas_hist)
    smax = maximum(salinity_meas_hist)
    spad = max(0.5, 0.05 * max(smax - smin, eps()))
    target_salinity = 30.0

    sal_cache = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
    # Pre-populate cache with the static first-frame data so animation uses the same static environment
    sal_cache[first_frame.file] = (lon_v_static, lat_v_static, salt_v_static)

    l = @layout [a b; c]

    anim = @animate for k in frame_ids
        lon_r, lat_r = lon_hist[k], lat_hist[k]
        
        tr_lon = lon_hist[1:k]
        tr_lat = lat_hist[1:k]

        fr = first_frame
        if !haskey(sal_cache, fr.file)
            sal_cache[fr.file] = load_surface_salinity(fr.file)
        end
        lon_v, lat_v, salt_v = sal_cache[fr.file]
    
        roi = (lon_v .>= xlims_map[1]) .& (lon_v .<= xlims_map[2]) .&
            (lat_v .>= ylims_map[1]) .& (lat_v .<= ylims_map[2])
        idx = findall(roi)
        
        if length(idx) > 12000
            step_val = cld(length(idx), 12000)
            idx = idx[1:step_val:end]
        end

        p1 = scatter(
            lon_v[idx], lat_v[idx];
            marker_z = salt_v[idx],
            m = (2.0, :square, stroke(0)),
            c = :viridis, clims = (0, 35),
            colorbar = true, colorbar_title = "Salinity (PSU)",
            xlabel = "Longitude", ylabel = "Latitude",
            title = "Boat + Salinity\n$(Dates.format(sim_times[k], dateformat"yyyy-mm-dd HH:MM")) UTC | Lap $(lap_hist[k])",
            aspect_ratio = :equal, xlims = xlims_map, ylims = ylims_map,
            framestyle = :box, grid = true, legend = :bottomleft
        )
        plot!(p1, transect_lon_shifted, transect_lat_shifted; lw = 2.5, lc = :white, label = "Planned Path")
        plot!(p1, tr_lon, tr_lat; lw = 2, lc = :dodgerblue, label = "Track")
        scatter!(p1, [lon_r], [lat_r]; ms = 6, mc = :red, label = "Boat")

        p2 = plot(
            lon_path, lat_path;
            line_z = clarity_history[:, k], 
            linewidth = 6, 
            c = :plasma, clims = (0, 1),
            colorbar = true, colorbar_title = "Clarity",
            xlabel = "Longitude", ylabel = "Latitude",
            title = "Clarity Map (Path Points)",
            aspect_ratio = :equal, xlims = xlims_map, ylims = ylims_map,
            legend = false, framestyle = :box, grid = true
        )
        
        scatter!(p2, [lon_r], [lat_r]; ms = 6, mc = :white, msw = 0, label = false)
        
        p3 = plot(
            sim_times[1:k], salinity_meas_hist[1:k];
            lw = 2.5, lc = :darkorange,
            xlabel = "Time", ylabel = "Measured Salinity (PSU)",
            title = "Measured Salinity vs Time",
            ylims = (smin - spad, smax + spad),
            xlims = (sim_times[1], sim_times[end]), 
            legend = :topright, framestyle = :box, grid = true,
            label = "Measured"
        )
        hline!(p3, [target_salinity]; lc = :red, ls = :dash, lw = 1.5, label = "Target")

        plot(p1, p2, p3; layout = l, size = (1400, 900))
    end

    mp4(anim, "static_env_speed_results_nobattery.mp4"; fps = 15)

    nsteps = length(sim_times)
    robot_point_idx = Vector{Int}(undef, nsteps)
    weight_at_robot = Vector{Float64}(undef, nsteps)

    for k in 1:nsteps
        d2 = (lon_path .- lon_hist[k]).^2 .+ (lat_path .- lat_hist[k]).^2
        idx = argmin(d2)
        robot_point_idx[k] = idx
        weight_at_robot[k] = weights_history[idx, k]
    end

    p = plot(
        sim_times, speed_hist;
        lw = 2.5, lc = :dodgerblue,
        xlabel = "Time",
        ylabel = "Speed (m/s)",
        title = "Speed vs Time with Robot-Point Weight Overlay",
        label = "Speed",
        framestyle = :box,
        grid = true,
        legend = :topright,
        ylims = (0, maximum(speed_hist)*1.1)
    )

    plot!(
        twinx(), sim_times, weight_at_robot;
        lw = 2.2, lc = :crimson, ls = :dash,
        ylabel = "Weight",
        label = "Weight at Robot Point"
    )
    savefig("static_env_speed_weight_overlay.png")


    final_stats = "✅ Success! Total Clarity: $(round(total_clarity, digits=2))"
    send_ntfy(final_stats, "Sim Finished", "default")
catch e
    error_log = "❌ Simulation Crashed: $(e)"
    send_ntfy(error_log, "Sim Failed", "high")
    rethrow(e)
end


# Run constant speed simulation for baseline comparison (using the same static environment)
function run_const_sim_static()
    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]

    # Precompute Local Flat-Earth Cartesian Coordinates for the Path
    R_earth = 6371000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    cos_mean_lat = cos(mean_lat_rad)

    X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
    Y_path = deg2rad.(lat_path) .* R_earth

    # Generate times
    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)

    # Initial clarity state (The physical transient starts at 0)
    clarity_state = zeros(Float64, npts)

    # History Arrays
    clarity_history    = zeros(Float64, npts, N_steps)
    weights_history    = zeros(Float64, npts, N_steps)
    salinity_meas_hist = zeros(Float64, N_steps)
    lon_hist           = zeros(Float64, N_steps)
    lat_hist           = zeros(Float64, N_steps)
    lap_hist           = zeros(Int, N_steps)
    speed_hist         = zeros(Float64, N_steps)

    # Buffers for the current lap's measurements
    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2000)
    sizehint!(lap_pos_buffer, 2000)

    s_robot_total = 0.0
    previous_lap = 1

    current_weights = fill(1.0, npts)

    # Use static environment already loaded: lon_v_static, lat_v_static, salt_v_static
    println("Starting Constant-Speed (Static Env) Simulation Loop...")

    @showprogress 1 "Running Const Sim..." for (t_idx, t) in enumerate(sim_times)

        # --- A. Determine Lap & Position ---
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m)

        # --- B. LAP TRANSITION LOGIC ---
        if current_lap > previous_lap
            # 1. Update weights based on collected lap data (Logged for comparison only)
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)

            # 2. Clear buffers for the next lap
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)

            previous_lap = current_lap
        end

        # --- C. Environmental Lookup & Buffering ---
        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        # Use static salinity arrays loaded once at startup
        lon_v, lat_v, salt_v = lon_v_static, lat_v_static, salt_v_static
        s_meas = sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)

        push!(lap_sal_buffer, s_meas)
        push!(lap_pos_buffer, s_mod)

        # =====================================================================
        # --- D. HIGH-RESOLUTION PHYSICS SWEEP ---
        # =====================================================================
        inner_dt = 0.25
        t_inner = 0.0
        path_length_true = s_vec[end]
        s_t = s_robot_total

        while t_inner < dt_step_sec
            dt_inc = min(inner_dt, dt_step_sec - t_inner)
            s_t_mod = mod(s_t, path_length_true)

            # --- FIX SPEED TO NOMINAL ---
            u_current = u_nominal

            X_t, Y_t = get_2d_position(s_t_mod, s_vec, X_path, Y_path)

            dx = X_path .- X_t
            dy = Y_path .- Y_t
            dist_sq = dx.^2 .+ dy.^2

            Sj_current = sensing_function.(dist_sq; S_0 = 1.0, sigma = sigma_val)
            clarity_state .= update_clarity_rk4.(clarity_state, Sj_current, dt_inc; alpha=alpha_val)

            s_t += u_current * dt_inc
            t_inner += dt_inc
        end

        s_robot_total = s_t

        # --- E. Log Everything ---
        clarity_history[:, t_idx]  = clarity_state
        weights_history[:, t_idx]  = current_weights
        salinity_meas_hist[t_idx]  = s_meas
        lon_hist[t_idx]            = lon_r
        lat_hist[t_idx]            = lat_r
        lap_hist[t_idx]            = current_lap
        speed_hist[t_idx]          = u_nominal
    end

    return clarity_history, speed_hist, lap_hist, lon_hist, lat_hist, salinity_meas_hist, weights_history
end

try
    println("Running constant-speed (static env) simulation for baseline comparison...")
    clarity_history_const, speed_hist_const, lap_hist_const, lon_hist_const, lat_hist_const, salinity_meas_hist_const, weights_history_const = run_const_sim_static()

    # Save results to JLD2 file for later analysis
    filename_const = Dates.format(now(), "yyyy-mm-dd_HHMMSS") * "_const_speed_static_sim_results.jld2"
    @save filename_const clarity_history_const speed_hist_const lap_hist_const lon_hist_const lat_hist_const salinity_meas_hist_const weights_history_const
    println("Constant-speed (static) simulation results saved to ", filename_const)

    # Compute and print summary statistic for comparison
    total_clarity_const = sum(weights_history_const .* clarity_history_const)
    println("Total Accumulated Clarity (Constant Speed, Static Env): ", round(total_clarity_const, digits=2))

    # Compare with optimized run (if available)
    if @isdefined total_clarity
        println("Optimized Total Clarity: ", round(total_clarity, digits=2))
        println("Difference (Opt - Const): ", round(total_clarity - total_clarity_const, digits=2))
    end

    # --------------------- Plotting & Analysis (Constant Speed) ---------------------
    npts_const = size(clarity_history_const, 1)
    lon_path_const, lat_path_const, s_vec_const = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts_const)
    sim_times_const = collect(frames[1].dt : Minute(1) : frames[end].dt)
    frame_times = [fr.dt for fr in frames]
    frame_idx_for_time_const(t) = argmin(abs.(Dates.value.(frame_times .- t)))

    total_laps_const = lap_hist_const[end]
    lap_clarity_const = zeros(Float64, total_laps_const)
    for lap in 1:total_laps_const
        idx_for_lap = findall(x -> x == lap, lap_hist_const)
        lap_clarity_const[lap] = sum(clarity_history_const[:, idx_for_lap])
        println("[Const] Lap $lap Total Clarity: ", round(lap_clarity_const[lap], digits=2))
    end

    target_salinity = 30.0
    salinity_deviation_const = salinity_meas_hist_const .- target_salinity

    println("[Const] Mean Salinity Deviation: ", round(mean(salinity_deviation_const), digits=2), " PSU")
    println("[Const] Std Dev Salinity Deviation: ", round(std(salinity_deviation_const), digits=2), " PSU")

    histogram(
        salinity_deviation_const,
        bins = 30,
        normalize = :probability,
        title = "Deviation from Target Salinity (30 PSU)\nConstant Speed (Static)",
        xlabel = "Deviation (PSU)",
        ylabel = "Frequency",
        label = "Measurements",
        color = :seagreen,
        linecolor = :white,
        grid = true,
        size = (600, 400)
    )
    vline!([0.0], lw=3, color=:red, label="Target (30 PSU)")
    savefig("const_static_env_speed_salinity_deviation_histogram.png")

    # Animation and maps
    nsteps_const = length(sim_times_const)
    @assert nsteps_const > 0 "No simulation history found for constant-speed run."

    frame_stride = max(1, fld(nsteps_const, 300))
    frame_ids = 1:frame_stride:nsteps_const

    lat_min, lat_max = 36.75, 37.5
    lon_pad, lat_pad = 0.02, 0.02
    xlims_map = (minimum(lon_path_const) - lon_pad, maximum(lon_path_const) + lon_pad)
    ylims_map = (lat_min - lat_pad, lat_max + lat_pad)

    smin = minimum(salinity_meas_hist_const)
    smax = maximum(salinity_meas_hist_const)
    spad = max(0.5, 0.05 * max(smax - smin, eps()))

    sal_cache_const = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
    sal_cache_const[first_frame.file] = (lon_v_static, lat_v_static, salt_v_static)

    l = @layout [a b; c]
    anim_const = @animate for k in frame_ids
        lon_r, lat_r = lon_hist_const[k], lat_hist_const[k]
        tr_lon = lon_hist_const[1:k]
        tr_lat = lat_hist_const[1:k]

        frc = first_frame
        lon_vc, lat_vc, salt_vc = sal_cache_const[frc.file]

        roi = (lon_vc .>= xlims_map[1]) .& (lon_vc .<= xlims_map[2]) .&
            (lat_vc .>= ylims_map[1]) .& (lat_vc .<= ylims_map[2])
        idx = findall(roi)
        if length(idx) > 12000
            step_val = cld(length(idx), 12000)
            idx = idx[1:step_val:end]
        end

        p1 = scatter(
            lon_vc[idx], lat_vc[idx];
            marker_z = salt_vc[idx],
            m = (2.0, :square, stroke(0)),
            c = :viridis, clims = (0, 35),
            colorbar = true, colorbar_title = "Salinity (PSU)",
            xlabel = "Longitude", ylabel = "Latitude",
            title = "Boat + Salinity\n$(Dates.format(sim_times_const[k], dateformat"yyyy-mm-dd HH:MM")) UTC | Lap $(lap_hist_const[k])",
            aspect_ratio = :equal, xlims = xlims_map, ylims = ylims_map,
            framestyle = :box, grid = true, legend = :bottomleft
        )
        plot!(p1, transect_lon_shifted, transect_lat_shifted; lw = 2.5, lc = :white, label = "Planned Path")
        plot!(p1, tr_lon, tr_lat; lw = 2, lc = :dodgerblue, label = "Track")
        scatter!(p1, [lon_r], [lat_r]; ms = 6, mc = :red, label = "Boat")

        p2 = plot(
            lon_path_const, lat_path_const;
            line_z = clarity_history_const[:, k],
            linewidth = 6,
            c = :plasma, clims = (0, 1),
            colorbar = true, colorbar_title = "Clarity",
            xlabel = "Longitude", ylabel = "Latitude",
            title = "Clarity Map (Constant Speed)",
            aspect_ratio = :equal, xlims = xlims_map, ylims = ylims_map,
            legend = false, framestyle = :box, grid = true
        )
        scatter!(p2, [lon_r], [lat_r]; ms = 6, mc = :white, msw = 0, label = false)

        p3 = plot(
            sim_times_const[1:k], salinity_meas_hist_const[1:k];
            lw = 2.5, lc = :darkorange,
            xlabel = "Time", ylabel = "Measured Salinity (PSU)",
            title = "Measured Salinity vs Time",
            ylims = (smin - spad, smax + spad),
            xlims = (sim_times_const[1], sim_times_const[end]),
            legend = :topright, framestyle = :box, grid = true,
            label = "Measured"
        )
        hline!(p3, [target_salinity]; lc = :red, ls = :dash, lw = 1.5, label = "Target")

        plot(p1, p2, p3; layout = l, size = (1400, 900))
    end

    mp4(anim_const, "const_static_env_speed_results_nobattery.mp4"; fps = 15)

    # Speed vs time + weight overlay
    nsteps = length(sim_times_const)
    robot_point_idx_const = Vector{Int}(undef, nsteps)
    weight_at_robot_const = Vector{Float64}(undef, nsteps)
    for k in 1:nsteps
        d2 = (lon_path_const .- lon_hist_const[k]).^2 .+ (lat_path_const .- lat_hist_const[k]).^2
        idx = argmin(d2)
        robot_point_idx_const[k] = idx
        weight_at_robot_const[k] = weights_history_const[idx, k]
    end

    p_const = plot(
        sim_times_const, speed_hist_const;
        lw = 2.5, lc = :dodgerblue,
        xlabel = "Time",
        ylabel = "Speed (m/s)",
        title = "Speed vs Time with Robot-Point Weight Overlay (Constant Speed)",
        label = "Speed",
        framestyle = :box,
        grid = true,
        legend = :topright,
        ylims = (0, maximum(speed_hist_const)*1.1)
    )

    plot!(
        twinx(), sim_times_const, weight_at_robot_const;
        lw = 2.2, lc = :crimson, ls = :dash,
        ylabel = "Weight",
        label = "Weight at Robot Point"
    )
    savefig("const_static_env_speed_weight_overlay.png")

    final_stats = "✅ Constant-speed (static) Success! Total Clarity: $(round(total_clarity_const, digits=2))"
    send_ntfy(final_stats, "Const Sim Finished", "default")
catch e
    error_log = "❌ Constant-speed (static) Simulation Crashed: $(e)"
    send_ntfy(error_log, "Const Sim Failed", "high")
    rethrow(e)
end

# Print statistics for comparison
if (@isdefined total_clarity) && (@isdefined total_clarity_const)
    println("\n=== SUMMARY COMPARISON ===")
    println("Optimized Total Clarity: ", round(total_clarity, digits=2))
    println("Constant-Speed Total Clarity: ", round(total_clarity_const, digits=2))
    println("Difference (Opt - Const): ", round(total_clarity - total_clarity_const, digits=2))
    println("Percent Improvement: ", round(100 * (total_clarity - total_clarity_const) / abs(total_clarity_const), digits=2), "%")
end