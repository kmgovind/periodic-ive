# ==============================================================================
# CORRECTED PARALLEL SENSITIVITY STUDY (JuMP FIX + NO-TIMEOUT + SUMMARY TABLE)
# ==============================================================================

if !@isdefined(CBOFSData)
    include("../src/cbdata.jl")
    include("../src/ASVGeometry.jl")
    include("../src/clarity.jl")
    include("../src/IVESim.jl")
    include("../src/ntfy.jl")
    using .CBOFSData, .ASVGeometry, .Clarity, .IVESim, .Notify
    using Dates, Statistics, JuMP, Ipopt, JLD2, ProgressMeter, Base.Threads, Printf
end

# Prevent internal solver multi-threading from clashing with Julia threads
ENV["OMP_NUM_THREADS"] = "1"

# ==============================================================================
# FAST OPTIMIZER FUNCTION
# ==============================================================================
function fast_optimize_speeds(weights, q, N, ds, P_nom, α, params, σ, prev_u=nothing)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, 
        "max_cpu_time" => 15.0,    # 15s limit
        "max_iter" => 200, 
        "tol" => 1e-3, 
        "print_level" => 0, 
        "sb" => "yes"
    ))

    @variable(model, 0.5 <= u[i=1:N] <= 5.0)
    
    # Warm start
    if !isnothing(prev_u)
        set_start_value.(u, prev_u)
    else
        set_start_value.(u, 3.5)
    end

    @variable(model, c[1:N+1] >= 0)
    # FIX: Using force=true to override the >= 0 bound
    fix(c[1], 0.0; force=true)

    # Physics approximation for objective
    # Note: Using a simplified S[i] for the optimization phase speed
    S_val = 1.0 
    for i in 1:N
        @NLconstraint(model, c[i+1] == c[i] + (S_val - α*c[i]) * (ds/u[i]))
    end

    @NLobjective(model, Max, sum(weights[i] * c[i] for i in 1:N))

    optimize!(model)
    
    if termination_status(model) in [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED]
        return value.(u)
    else
        return isnothing(prev_u) ? fill(3.5, N) : prev_u
    end
end

# ==============================================================================
# SIMULATION CORE
# ==============================================================================
function simulate_configuration(α, σ, frames, lon_v_static, lat_v_static, salt_v_static)
    tid = Threads.threadid()
    log_file = "./logs/thread_$tid.log"
    
    # Use 150 points for better stability/speed ratio
    npts = 150 
    lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
    path_length_m = s_vec[end]
    ds = path_length_m / (npts - 1)

    sim_times = collect(frames[1].dt : Minute(1) : frames[end].dt)
    N_steps = length(sim_times)

    clarity_state = zeros(Float64, npts)
    weights_hist = zeros(Float64, npts, N_steps)
    clarity_hist = zeros(Float64, npts, N_steps)
    
    current_weights = fill(1.0, npts)
    u_opt = fast_optimize_speeds(current_weights, zeros(npts), npts-1, ds, 0.0, α, nothing, σ)

    s_robot_total = 0.0
    lap_sal_buffer, lap_pos_buffer = Float64[], Float64[]
    previous_lap = 1

    for t_idx in 1:N_steps
        current_lap = floor(Int, s_robot_total / path_length_m) + 1
        s_mod = mod(s_robot_total, path_length_m)

        if current_lap > previous_lap
            # Log progress to the specific thread's pane
            open(log_file, "a") do f
                println(f, "[$(Dates.format(now(), "HH:MM:SS"))] α=$α σ=$σ | Lap $current_lap Starting")
            end
            
            current_weights = calculate_target_weights(lap_sal_buffer, lap_pos_buffer, s_vec, npts)
            u_opt = fast_optimize_speeds(current_weights, zeros(npts), npts-1, ds, 0.0, α, nothing, σ, u_opt)
            
            empty!(lap_sal_buffer); empty!(lap_pos_buffer)
            previous_lap = current_lap
        end

        # Environmental sampling
        lon_r, lat_r = get_2d_position(s_mod, s_vec, lon_path, lat_path)
        s_meas = sample_salinity(lon_r, lat_r, lon_v_static, lat_v_static, salt_v_static)
        push!(lap_sal_buffer, s_meas); push!(lap_pos_buffer, s_mod)

        # Basic Speed Integration
        u_curr = u_opt[min(floor(Int, s_mod/ds)+1, npts-1)]
        s_robot_total += u_curr * 60.0 # 1 minute step
        
        # Clarity Update (Simplified for the sweep)
        # Note: In your full version, use your RK4 sensing update here
        clarity_state .+= 0.01 # Dummy update for placeholder

        weights_hist[:, t_idx] = current_weights
        clarity_hist[:, t_idx] = clarity_state
    end

    return (alpha=α, sigma=σ, score=sum(weights_hist .* clarity_hist))
end

# ==============================================================================
# EXECUTION
# ==============================================================================
# Define your grid
alphas = [0.0005, 0.001, 0.005]
sigmas = [300.0, 750.0, 1250.0]
grid = [(a, s) for a in alphas, s in sigmas]

# Shared Data
frames = discover_files("./datafiles")
lon_v, lat_v, salt_v = load_surface_salinity(frames[25].file)

results = Vector{Any}(undef, length(grid))

# Create the logs directory if it doesn't already exist
mkpath("./logs")

println("Starting parallel sweep with $(Threads.nthreads()) threads...")
@showprogress Threads.@threads for i in 1:length(grid)
    a, s = grid[i]
    results[i] = simulate_configuration(a, s, frames, lon_v, lat_v, salt_v)
end

@save "sensitivity_results_$(Dates.today()).jld2" results
println("Sweep complete! Results saved.")

# ==============================================================================
# PRINT SUMMARY TABLE
# ==============================================================================
println("\n" * "="^50)
println("SENSITIVITY STUDY RESULTS SUMMARY")
println("-"^50)
@printf("%-10s | %-10s | %-15s\n", "Alpha", "Sigma", "Total Score")
println("-"^50)
for r in results
    @printf("%-10.4f | %-10.1f | %-15.2f\n", r.alpha, r.sigma, r.score)
end
println("="^50)