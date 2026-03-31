# ==============================================================================
# VARIABLE & CONSTANT SPEED SIMULATION ENGINE (DYNAMIC ENVIRONMENT)
# ==============================================================================

if !isdefined(@__MODULE__, :CBOFSData)
    println("Loading modules and packages...")
    include("../src/cbdata.jl")
    include("../src/ASVGeometry.jl")
    include("../src/clarity.jl")
    include("../src/IVESim.jl")
    include("../src/ntfy.jl")
    
    using .CBOFSData, .ASVGeometry, .Clarity, .IVESim, .Notify
    using Dates, Statistics, JuMP, Ipopt, JLD2, ProgressMeter, Base.Threads, Printf
end

# --- Global Simulation Parameters ---
u_nominal = 1.75 # m/s
dt_step_sec = 60.0 # seconds
sigma_val = 1500.0 # sensing footprint in meters (Adjusted for large spatial runs!)
alpha_val = 0.001 # clarity decay rate

# ==============================================================================
# 1. OPTIMIZED VARIABLE-SPEED SIMULATION
# ==============================================================================
function run_dynamic_opt_sim(frames)
    vehicle_params = VehicleParams(u_nominal, 24*3600.0)

    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]

    # Precompute Local Flat-Earth Cartesian Coordinates
    R_earth = 6371000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    cos_mean_lat = cos(mean_lat_rad)
    X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
    Y_path = deg2rad.(lat_path) .* R_earth

    P_in_W = vehicle_params.kh + vehicle_params.km * (u_nominal^3)
    N_segments = npts - 1
    ds = path_length_m / N_segments

    # Generate times 
    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)
    
    # Frame lookup helpers
    frame_times = [fr.dt for fr in frames]
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(frame_times .- t)))
    current_frame_idx = -1
    local lon_v, lat_v, salt_v

    clarity_state = zeros(Float64, npts)
    clarity_history    = zeros(Float64, npts, N_steps)
    weights_history    = zeros(Float64, npts, N_steps)
    salinity_meas_hist = zeros(Float64, N_steps)
    lon_hist           = zeros(Float64, N_steps)
    lat_hist           = zeros(Float64, N_steps)
    lap_hist           = zeros(Int, N_steps)
    speed_hist         = zeros(Float64, N_steps)

    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2000) 
    sizehint!(lap_pos_buffer, 2000)

    s_robot_total = 0.0 
    previous_lap = 1 

    # --- LAP 1 OFFLINE OPTIMIZATION ---
    current_weights = fill(1.0, npts)
    q_initial = fill(0.0, npts) 
    u_initial = fill(u_nominal, N_segments) 
    
    println("Running Lap 1 Offline Optimization...")
    u_opt_segments = optimize_lap_speeds_full_spatiotemporal(
        current_weights, q_initial, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val
    )

    # --- MAIN SIMULATION LOOP ---
    println("Starting Optimized Variable-Speed Simulation (Dynamic Env)...")
    @showprogress 1 "Opt Sim..." for (t_idx, t) in enumerate(sim_times)
        
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m) 
        
        if current_lap > previous_lap
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)
            u_initial = copy(u_opt_segments)
            u_opt_segments = optimize_lap_speeds_full_spatiotemporal_warm(
                current_weights, q_initial, u_initial, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val
            )
            empty!(lap_sal_buffer)
            empty!(lap_pos_buffer)
            previous_lap = current_lap
        end
        
        # --- Environmental Lookup (Smart Cache) ---
        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        
        f_idx = frame_idx_for_time(t)
        if f_idx != current_frame_idx
            lon_v, lat_v, salt_v = load_surface_salinity(frames[f_idx].file)
            current_frame_idx = f_idx
        end
        
        s_meas = sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)
        push!(lap_sal_buffer, s_meas)
        push!(lap_pos_buffer, s_mod)

        # --- HIGH-RESOLUTION PHYSICS SWEEP ---
        inner_dt = 0.25 
        t_inner = 0.0
        s_t = s_robot_total
        u_current_display = 0.0
        
        while t_inner < dt_step_sec
            dt_inc = min(inner_dt, dt_step_sec - t_inner)
            s_t_mod = mod(s_t, path_length_m)
            
            seg_idx = min(floor(Int, s_t_mod / ds) + 1, N_segments)
            u_current = u_opt_segments[seg_idx]
            u_current_display = u_current
            
            X_t, Y_t = get_2d_position(s_t_mod, s_vec, X_path, Y_path)
            dx = X_path .- X_t
            dy = Y_path .- Y_t
            dist_sq = dx.^2 .+ dy.^2
            
            Sj_current = sensing_function.(dist_sq; S_0 = 0.005, sigma = sigma_val)
            clarity_state .= update_clarity_rk4.(clarity_state, Sj_current, dt_inc; alpha=alpha_val)
            
            s_t += u_current * dt_inc
            t_inner += dt_inc
        end
        
        s_robot_total = s_t
        
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


# ==============================================================================
# 2. CONSTANT SPEED BASELINE SIMULATION
# ==============================================================================
function run_dynamic_const_sim(frames)
    npts = 200
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]

    R_earth = 6371000.0
    mean_lat_rad = deg2rad(mean(lat_path))
    cos_mean_lat = cos(mean_lat_rad)
    X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
    Y_path = deg2rad.(lat_path) .* R_earth

    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)

    frame_times = [fr.dt for fr in frames]
    frame_idx_for_time(t) = argmin(abs.(Dates.value.(frame_times .- t)))
    current_frame_idx = -1
    local lon_v, lat_v, salt_v

    clarity_state = zeros(Float64, npts)
    clarity_history    = zeros(Float64, npts, N_steps)
    weights_history    = zeros(Float64, npts, N_steps)
    salinity_meas_hist = zeros(Float64, N_steps)
    lon_hist           = zeros(Float64, N_steps)
    lat_hist           = zeros(Float64, N_steps)
    lap_hist           = zeros(Int, N_steps)
    speed_hist         = zeros(Float64, N_steps)

    lap_sal_buffer = Float64[]
    lap_pos_buffer = Float64[]
    sizehint!(lap_sal_buffer, 2000)
    sizehint!(lap_pos_buffer, 2000)

    s_robot_total = 0.0
    previous_lap = 1
    current_weights = fill(1.0, npts)

    println("Starting Constant-Speed Baseline (Dynamic Env)...")
    @showprogress 1 "Const Sim..." for (t_idx, t) in enumerate(sim_times)

        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m)

        if current_lap > previous_lap
            # Update weights purely for logging and fair evaluation (doesn't affect speed)
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)
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

        # --- HIGH-RESOLUTION PHYSICS SWEEP ---
        inner_dt = 0.25
        t_inner = 0.0
        s_t = s_robot_total

        while t_inner < dt_step_sec
            dt_inc = min(inner_dt, dt_step_sec - t_inner)
            s_t_mod = mod(s_t, path_length_m)
            u_current = u_nominal # CONSTANT SPEED ENFORCED HERE

            X_t, Y_t = get_2d_position(s_t_mod, s_vec, X_path, Y_path)
            dx = X_path .- X_t
            dy = Y_path .- Y_t
            dist_sq = dx.^2 .+ dy.^2

            Sj_current = sensing_function.(dist_sq; S_0 = 0.005, sigma = sigma_val)
            clarity_state .= update_clarity_rk4.(clarity_state, Sj_current, dt_inc; alpha=alpha_val)

            s_t += u_current * dt_inc
            t_inner += dt_inc
        end

        s_robot_total = s_t

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

# ==============================================================================
# 3. MAIN EXECUTION BLOCK
# ==============================================================================
frames = discover_files("./datafiles")

# --- Run & Save Optimized Version ---
try
    clarity_opt, speed_opt, lap_opt, lon_opt, lat_opt, sal_opt, weights_opt = run_dynamic_opt_sim(frames)

    filename_opt = Dates.format(now(), "yyyy-mm-dd_HHMMSS") * "_dynamic_opt_sim.jld2"
    @save filename_opt clarity_history=clarity_opt speed_hist=speed_opt lap_hist=lap_opt lon_hist=lon_opt lat_hist=lat_opt salinity_meas_hist=sal_opt weights_history=weights_opt
    println("✅ Optimized dynamic results saved to ", filename_opt)
    
    global total_clarity_opt = sum(weights_opt .* clarity_opt)
catch e
    send_ntfy("❌ Optimized Dynamic Sim Crashed: $(e)", "Sim Failed", "high")
    rethrow(e)
end

# --- Run & Save Constant Baseline Version ---
try
    clarity_const, speed_const, lap_const, lon_const, lat_const, sal_const, weights_const = run_dynamic_const_sim(frames)

    filename_const = Dates.format(now(), "yyyy-mm-dd_HHMMSS") * "_dynamic_const_sim.jld2"
    @save filename_const clarity_history=clarity_const speed_hist=speed_const lap_hist=lap_const lon_hist=lon_const lat_hist=lat_const salinity_meas_hist=sal_const weights_history=weights_const
    println("✅ Constant-speed dynamic results saved to ", filename_const)
    
    global total_clarity_const = sum(weights_const .* clarity_const)
catch e
    send_ntfy("❌ Constant Dynamic Sim Crashed: $(e)", "Sim Failed", "high")
    rethrow(e)
end

# --- Print Final Comparison Summary ---
if (@isdefined total_clarity_opt) && (@isdefined total_clarity_const)
    println("\n" * "="^40)
    println("📊 DYNAMIC SIMULATION SUMMARY")
    println("="^40)
    println("Optimized Total Clarity:     ", round(total_clarity_opt, digits=2))
    println("Constant-Speed Total Clarity:", round(total_clarity_const, digits=2))
    println("Difference (Opt - Const):    ", round(total_clarity_opt - total_clarity_const, digits=2))
    println("Percent Improvement:         ", round(100 * (total_clarity_opt - total_clarity_const) / abs(total_clarity_const), digits=2), "%")
    println("="^40)
    
    send_ntfy("✅ Dynamic sims finished! Opt: $(round(total_clarity_opt, digits=1)) | Const: $(round(total_clarity_const, digits=1))", "Sims Complete", "default")
end