# Only include and load packages if they haven't been loaded yet in this session
if !@isdefined(CBOFSData)
    println("Loading modules and packages...")
    include("../src/cbdata.jl")
    include("../src/ASVGeometry.jl")
    include("../src/clarity.jl")

    using .CBOFSData, .ASVGeometry, .Clarity, .Notify
    using Dates, Statistics, Plots, JuMP, Ipopt, ProgressMeter
    
    # Optional: If you are actively editing cbdata.jl, ASVGeometry.jl, etc. 
    # you can use Revise to automatically track changes without re-including!
    using Revise
end


################### Transect Definition #######################################################################################
### Load salinity data
frames = discover_files("./datafiles")

fr = frames[1]
lon_valid, lat_valid, salt_valid = load_surface_salinity(fr.file)
roi = (lat_valid .>= lat_min) .& (lat_valid .<= lat_max)

try

    ###############################  Simulation   ###########################################
    # Initialize ASV geometry and simulate path discretization
    u_nominal = 1.75  # m/s
    vehicle_params = VehicleParams(u_nominal, 24*3600.0)  # 24 hours at 1.75 m/s nominal speed

    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]

    # Precompute Local Flat-Earth Cartesian Coordinates for the Path
    const R_earth = 6371000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    cos_mean_lat = cos(mean_lat_rad)

    # X and Y coordinates in meters
    X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
    Y_path = deg2rad.(lat_path) .* R_earth


    # Simulation parameters
    dt_step_sec = 60.0              # 1 minute per step
    sigma_val = 3000.0              # Sensing footprint [m]
    alpha_val = 0.05                # Clarity decay

    # NEW: Calculate Power Budget based on nominal baseline speed
    power_nominal = vehicle_params.kh + vehicle_params.km * (u_nominal^3)
    P_in_W = power_nominal          # Continuous Solar Power Generation [Watts]
    # P_in_W = 750.0

    # Path Discretization for JuMP
    N_segments = npts - 1               # Number of path segments
    ds = path_length_m / N_segments     # Length of each segment [m]

    # Generate times (assuming frames are 1 hour apart, step every minute)
    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)

    # Initial clarity state
    clarity_state = zeros(Float64, npts)

    # History Arrays
    clarity_history    = zeros(Float64, npts, N_steps)
    weights_history    = zeros(Float64, npts, N_steps)
    salinity_meas_hist = zeros(Float64, N_steps)
    lon_hist           = zeros(Float64, N_steps)
    lat_hist           = zeros(Float64, N_steps)
    battery_hist       = zeros(Float64, N_steps) 
    lap_hist           = zeros(Int, N_steps)
    speed_hist         = zeros(Float64, N_steps)

    # Buffers for the current lap's measurements
    # Pre-allocate memory for ~1.5 days worth of minutes just to be safe
    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2000) 
    sizehint!(lap_pos_buffer, 2000)

    s_robot_total = 0.0 
    previous_lap = 1                       

    # NEW: Battery Capacity Setup (e.g., capable of holding 24 hours of solar power)
    # battery_capacity_Wh = P_in_W * 24.0 
    battery_capacity_Wh = 6500.0
    current_battery = battery_capacity_Wh 

    frame_times = [fr.dt for fr in frames]
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(frame_times .- t)))


    # Run Initial Optimization for Lap 1 (Uniform weights)
    # low_w  = 0.2
    # high_w = 1.0

    # s_start = (4/8) * path_length_m
    # s_end   = (6/8) * path_length_m

    # current_weights = fill(low_w, npts)
    # current_weights[(s_vec .>= s_start) .& (s_vec .<= s_end)] .= high_w
    current_weights = fill(1.0, npts)
    opt_time_start = time()
    u_opt_segments = optimize_lap_speeds_full_spatiotemporal(current_weights, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val)
    opt_time_elapsed = time() - opt_time_start
    println("Initial optimization completed in $(round(opt_time_elapsed, digits=2)) seconds")


    # ==============================================================================
    # 4. MAIN SIMULATION LOOP
    # ==============================================================================
    println("Starting Optimized Variable-Speed Simulation Loop...")

    @showprogress 1 "Running Sim..." for (t_idx, t) in enumerate(sim_times)

        global s_robot_total, current_battery, previous_lap, current_weights, u_opt_segments
        
        # --- A. Determine Lap & Position ---
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m) 
        
        # --- B. LAP TRANSITION LOGIC ---
        if current_lap > previous_lap
            # 1. Calculate new weights from the lap that just finished
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)
            
            # 2. Run Optimization for the new lap 
            u_opt_segments = optimize_lap_speeds_full_spatiotemporal(current_weights, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val)
            
            # 3. Clear buffers for the new lap
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            
            previous_lap = current_lap
        end
        
        # --- C. Environmental Lookup & Buffering ---
        # Sample salinity once at the start of the 60-second frame
        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        frame_idx = frame_idx_for_time(t)
        lon_v, lat_v, salt_v = load_surface_salinity(frames[frame_idx].file) 
        s_meas = sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)
        
        push!(lap_sal_buffer, s_meas)
        push!(lap_pos_buffer, s_mod)

        # =====================================================================
        # --- D. COMBINED HIGH-RESOLUTION PHYSICS SWEEP ---
        # (Absorbs previous Sections C, E, F, and H)
        # =====================================================================
        inner_dt = 0.25 # Finer time steps for accurate ODE integration & physics
        t_inner = 0.0
        R_earth = 6371000.0
        
        # Use the true array length to prevent wrap-around coordinate teleports
        path_length_true = s_vec[end]
        
        # Local trackers for the continuous 60-second sweep
        s_t = s_robot_total
        bat_t = current_battery
        u_current_display = 0.0
        
        while t_inner < dt_step_sec
            dt_inc = min(inner_dt, dt_step_sec - t_inner)
            
            # 1. Exact position modulo true path length
            s_t_mod = mod(s_t, path_length_true)
            
            # 2. DYNAMIC SPEED LOOKUP (Fixes the battery desync)
            # Evaluated every 0.25s, so the robot changes speed the millisecond it crosses a segment!
            seg_idx = min(floor(Int, s_t_mod / ds) + 1, N_segments)
            u_opt = u_opt_segments[seg_idx]
            
            # 3. HIGH-RES BATTERY PHYSICS (Constant Power, Unbounded Capacity)
            energy_gained_Wh = P_in_W * (dt_inc / 3600.0)
            
            power_req_opt = vehicle_params.kh + vehicle_params.km * (u_opt^3)
            energy_req_opt_Wh = power_req_opt * (dt_inc / 3600.0)
            
            # Check if we have enough battery to maintain optimal speed
            if bat_t + energy_gained_Wh >= energy_req_opt_Wh
                u_current = u_opt
                bat_t = bat_t + energy_gained_Wh - energy_req_opt_Wh
            else
                # Stall out if budget is blown (should rarely happen with the equality constraint)
                u_current = 0.0
                energy_idle_Wh = vehicle_params.kh * (dt_inc / 3600.0)
                if bat_t + energy_gained_Wh >= energy_idle_Wh
                    bat_t = bat_t + energy_gained_Wh - energy_idle_Wh
                else
                    bat_t = 0.0 # Floor it at zero so we don't get negative energy
                end
            end
            
            # REMOVED: bat_t = min(battery_capacity_Wh, bat_t)
            # The battery can now buffer as much energy as it needs to without artificially capping
            u_current_display = u_current
            
            # 4. CLARITY SWEEP (Optimized Cartesian Distance)
            # Interpolate directly in meter-space to avoid all trigonometry
            X_t, Y_t = get_2d_position(s_t_mod, s_vec, X_path, Y_path)
            
            for j in 1:npts
                # Simple Euclidean distance in pre-computed meters
                dx = X_path[j] - X_t
                dy = Y_path[j] - Y_t
                dist_meters = sqrt(dx^2 + dy^2)
                
                # (Optional minor speedup: if your sensing_function can be modified to 
                # accept distance-squared, you can delete the sqrt() above entirely!)
                Sj_current = sensing_function(dist_meters; S_0 = 1.0, sigma = sigma_val)
                clarity_state[j] = update_clarity_rk4(clarity_state[j], Sj_current, dt_inc; alpha=alpha_val)
            end
            
            # 5. MOVE ROBOT 
            s_t += u_current * dt_inc
            t_inner += dt_inc
        end
        
        # Sync local variables back to the main loop state at the end of the 60s
        s_robot_total = s_t
        current_battery = bat_t
        u_current = u_current_display # Pass out the final speed for logging
        
        # --- E. Log Everything ---
        clarity_history[:, t_idx]  = clarity_state
        weights_history[:, t_idx]  = current_weights
        salinity_meas_hist[t_idx]  = s_meas
        lon_hist[t_idx]            = lon_r
        lat_hist[t_idx]            = lat_r
        battery_hist[t_idx]        = current_battery 
        lap_hist[t_idx]            = current_lap
        speed_hist[t_idx]          = u_current
    end


    ###################################################################################################



    ################################### Plotting & Analysis ########################################
    # --- Analytics 1: Total Clarity ---
    # Total clarity accumulated across all points and all times
    total_clarity = sum(weights_history .* clarity_history)
    println("Total Accumulated Clarity (All Time): ", round(total_clarity, digits=2))

    # --- Analytics 2: Per-Lap Clarity ---
    total_laps = lap_hist[end]
    lap_clarity = zeros(Float64, total_laps)

    for lap in 1:total_laps
        # Find all time indices belonging to this lap
        idx_for_lap = findall(x -> x == lap, lap_hist)
        # Sum the clarity matrix just for those time columns
        lap_clarity[lap] = sum(clarity_history[:, idx_for_lap])
        println("Lap $lap Total Clarity: ", round(lap_clarity[lap], digits=2))
    end

    # --- Analytics 3: Salinity Deviation Histogram ---
    target_salinity = 30.0
    # Calculate the deviation (Measurement - Target)
    salinity_deviation = salinity_meas_hist .- target_salinity

    println("Mean Salinity Deviation: ", round(mean(salinity_deviation), digits=2), " PSU")
    println("Standard Deviation of Salinity Deviation: ", round(std(salinity_deviation), digits=2), " PSU")

    # Plot the Histogram
    histogram(
        salinity_deviation,
        bins = 30,
        normalize = :probability,
        title = "Deviation from Target Salinity (30 PSU)\nVariable Speed",
        xlabel = "Deviation (PSU)",
        ylabel = "Frequency",
        label = "Measurements",
        color = :steelblue,
        linecolor = :white,
        grid = true,
        size = (600, 400)
    )
    # Add a vertical line at 0 (Target)
    vline!([0.0], lw=3, color=:red, label="Target (30 PSU)")
    savefig("variable_speed_salinity_deviation_histogram.png")


    nsteps = length(sim_times)
    @assert nsteps > 0 "No simulation history found. Run the simulation cell first."

    # Keep video length manageable (e.g., skip frames if the simulation is days long)
    frame_stride = max(1, fld(nsteps, 300)) # Targets around ~300 frames total
    frame_ids = 1:frame_stride:nsteps

    lat_min, lat_max = 36.75, 37.5
    lon_pad, lat_pad = 0.02, 0.02
    xlims_map = (minimum(lon_path) - lon_pad, maximum(lon_path) + lon_pad)
    ylims_map = (lat_min - lat_pad, lat_max + lat_pad)
    bmax = max(maximum(battery_hist), 1.0)

    smin = minimum(salinity_meas_hist)
    smax = maximum(salinity_meas_hist)
    spad = max(0.5, 0.05 * max(smax - smin, eps()))
    target_salinity = 30.0

    # Simple cache so we don't reload the NetCDF file 60 times for the same hour
    sal_cache = Dict{String, Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()

    # --- 3. Animation Loop ---

    anim = @animate for k in frame_ids
        lon_r, lat_r = lon_hist[k], lat_hist[k]
        
        # Track history up to current frame
        tr_lon = lon_hist[1:k]
        tr_lat = lat_hist[1:k]

        # Salinity background for current time
        fr = frames[frame_idx_for_time(sim_times[k])]
        if !haskey(sal_cache, fr.file)
            sal_cache[fr.file] = load_surface_salinity(fr.file)
        end
        lon_v, lat_v, salt_v = sal_cache[fr.file]
        
        # Filter background to ROI to speed up plotting
        roi = (lon_v .>= xlims_map[1]) .& (lon_v .<= xlims_map[2]) .&
            (lat_v .>= ylims_map[1]) .& (lat_v .<= ylims_map[2])
        idx = findall(roi)
        
        # Subsample background if it's too dense
        if length(idx) > 12000
            step_val = cld(length(idx), 12000)
            idx = idx[1:step_val:end]
        end

        # --- Panel 1: Boat position with salinity background ---
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

        # --- Panel 2: Clarity (Using continuous line instead of scatter) ---
        p2 = plot(
            lon_path, lat_path;
            line_z = clarity_history[:, k],   # Changed from marker_z
            linewidth = 6,                    # Set how thick you want the path
            c = :plasma, clims = (0, 1),
            colorbar = true, colorbar_title = "Clarity",
            xlabel = "Longitude", ylabel = "Latitude",
            title = "Clarity Map (Path Points)",
            aspect_ratio = :equal, xlims = xlims_map, ylims = ylims_map,
            legend = false, framestyle = :box, grid = true
        )
        
        # We still use scatter! for the robot itself, since it's just one point
        scatter!(p2, [lon_r], [lat_r]; ms = 6, mc = :white, msw = 0, label = false)
        
        # --- Panel 3: Battery vs time ---
        p3 = plot(
            sim_times[1:k], battery_hist[1:k];
            lw = 2.5, lc = :seagreen4,
            xlabel = "Time", ylabel = "Battery (Wh)",
            title = "Battery vs Time",
            ylims = (0, 1.05 * bmax),
            xlims = (sim_times[1], sim_times[end]), # Fix X axis so it doesn't jump around
            legend = false, framestyle = :box, grid = true
        )

        # --- Panel 4: Measured salinity vs time ---
        p4 = plot(
            sim_times[1:k], salinity_meas_hist[1:k];
            lw = 2.5, lc = :darkorange,
            xlabel = "Time", ylabel = "Measured Salinity (PSU)",
            title = "Measured Salinity vs Time",
            ylims = (smin - spad, smax + spad),
            xlims = (sim_times[1], sim_times[end]), # Fix X axis
            legend = :topright, framestyle = :box, grid = true,
            label = "Measured"
        )
        hline!(p4, [target_salinity]; lc = :red, ls = :dash, lw = 1.5, label = "Target")

        # Combine into a 2x2 grid
        plot(p1, p2, p3, p4; layout = (2, 2), size = (1500, 900))
    end

    mp4(anim, "variable_speed_results.mp4"; fps = 15)



    # Weight of the path point closest to the robot at each time step
    nsteps = length(sim_times)
    robot_point_idx = Vector{Int}(undef, nsteps)
    weight_at_robot = Vector{Float64}(undef, nsteps)

    for k in 1:nsteps
        d2 = (lon_path .- lon_hist[k]).^2 .+ (lat_path .- lat_hist[k]).^2
        idx = argmin(d2)
        robot_point_idx[k] = idx
        weight_at_robot[k] = weights_history[idx, k]
    end

    # Plot speed and overlaid weight (secondary y-axis)
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
        ylims = (0, 2.5)
    )

    plot!(
        twinx(), sim_times, weight_at_robot;
        lw = 2.2, lc = :crimson, ls = :dash,
        ylabel = "Weight",
        label = "Weight at Robot Point"
    )
    savefig("speed_weight_overlay.png")


    final_stats = "Success! \nTotal Clarity: $(round(total_clarity, digits=2))\nLaps Completed: $(lap_hist[end])"
    send_ntfy(final_stats, "Sim Finished ✅", "default")
catch e
    error_log = "Simulation Crashed. \nError: $(replace(string(e), "\"" => "'"))"
    send_ntfy(error_log, "Sim Failed ❌", "high")
    rethrow(e)
end