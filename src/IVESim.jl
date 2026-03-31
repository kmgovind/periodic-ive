module IVESim

# include("clarity.jl")

using Ipopt, JuMP, ..Clarity

export optimize_lap_speeds_full_spatiotemporal, optimize_lap_speeds_full_spatiotemporal_warm, get_2d_position

function optimize_lap_speeds_full_spatiotemporal_warm(
    weights, q_initial, u_initial, N, ds, P_in_W_avg, alpha_decay, vehicle_params, sigma_val
)

    model = Model(Ipopt.Optimizer)
    
    # set_silent(model)

    set_attribute(model, "tol", 1e-4)
    set_attribute(model, "acceptable_tol", 1e-3)
    set_attribute(model, "max_iter", 1000) 
    set_attribute(model, "max_cpu_time", 600.0)
    set_attribute(model, "mumps_mem_percent", 10000)
    set_optimizer_attribute(model, "warm_start_init_point", "yes")
    # set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
    
    S_matrix = zeros(N, N)
    path_length = N * ds
    for j in 1:N
        for k in 1:N
            dist_direct = abs((j-1)*ds - (k-1)*ds)
            dist_wrap = path_length - dist_direct
            shortest_dist = min(dist_direct, dist_wrap)
            val = sensing_function(shortest_dist^2; S_0=0.005, sigma=sigma_val)
            S_matrix[j, k] = val > 1e-5 ? val : 0.0
        end
    end
    
    # --- PHYSICAL TIME BOUNDS (Change of Variables!) ---
    safe_u_min = max(0.5, vehicle_params.u_min) 
    safe_u_max = 5.0 
    
    dt_min = ds / safe_u_max
    dt_max = ds / safe_u_min
    
    # Optimize OVER TIME instead of speed to eliminate 1/u nonlinearity!
    @variable(model, dt_min <= dt_var[1:N] <= dt_max, start = ds/1.75)
    
    # Speed is now just an expression (used only for energy)
    @NLexpression(model, u_expr[k=1:N], ds / dt_var[k])
    
    # T_lap is now a purely linear sum
    @expression(model, T_lap, sum(dt_var[k] for k in 1:N))
    
    @variable(model, 0 <= q[1:N, 1:N+1] <= 1.0, start=0.5)
    @variable(model, node_clarity[1:N] >= 0)
    
    # 1. Calculate the true Time-Integral of clarity for each node
    for j in 1:N
        @NLconstraint(model, 
            node_clarity[j] == sum( q[j, k+1] * dt_var[k] for k in 1:N )
        )
    end
    
    # 2. THE EPIGRAPH VARIABLE
    @variable(model, avg_clarity >= 0)
    
    # 3. THE EPIGRAPH CONSTRAINT
    @NLconstraint(model, 
        avg_clarity * T_lap <= sum(weights[j] * node_clarity[j] for j in 1:N)
    )
    
    # 4. PURELY LINEAR OBJECTIVE
    @objective(model, Max, avg_clarity)

    # Energy Constraint
    scale_factor = 1e6
    @NLconstraint(model, 
        sum( (vehicle_params.kh * dt_var[k] + vehicle_params.km * (u_expr[k]^2) * ds) / scale_factor for k in 1:N) <= (P_in_W_avg * T_lap) / scale_factor
    )
    
    # --- EXACT NONLINEAR DYNAMICS (Implicit Euler) ---
    for k in 1:N
        k_next = (k == N) ? 1 : k + 1 
        for j in 1:N
            S_k = S_matrix[j, k]
            S_k_next = S_matrix[j, k_next]
            
            if S_k == 0.0 && S_k_next == 0.0
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] - dt_var[k] * alpha_decay * q[j, k+1]^2
                )
            else
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] + dt_var[k] * (
                        S_k_next * (1 - q[j, k+1])^2 - alpha_decay * q[j, k+1]^2
                    )
                )
            end
        end
    end

    # Slack Periodic Boundaries
    @constraint(model, [j=1:N], q[j, 1] - q[j, N+1] <= 0.01)
    @constraint(model, [j=1:N], q[j, N+1] - q[j, 1] <= 0.01)

    # Warm start for times based on initial speed guess
    for i in 1:N
        set_start_value(dt_var[i], ds / u_initial[i])
    end
    
    optimize!(model)
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]
        println("Success! Speeds optimized.")
        return [ds / v for v in value.(dt_var)] 
    else
        println("Warning: Optimizer failed. Status: ", termination_status(model))
        return fill(1.75, N) 
    end
end

function optimize_lap_speeds_full_spatiotemporal(
    weights, q_initial, N, ds, P_in_W_avg, alpha_decay, vehicle_params, sigma_val
)
    model = Model(Ipopt.Optimizer)
    
    set_silent(model)

    set_attribute(model, "tol", 1e-4)
    set_attribute(model, "acceptable_tol", 1e-3)
    set_attribute(model, "max_iter", 1000) 
    set_attribute(model, "max_cpu_time", 600.0)
    set_attribute(model, "mumps_mem_percent", 10000)
    set_optimizer_attribute(model, "warm_start_init_point", "yes")
    
    S_matrix = zeros(N, N)
    path_length = N * ds
    for j in 1:N
        for k in 1:N
            dist_direct = abs((j-1)*ds - (k-1)*ds)
            dist_wrap = path_length - dist_direct
            shortest_dist = min(dist_direct, dist_wrap)
            val = sensing_function(shortest_dist^2; S_0=0.005, sigma=sigma_val)
            S_matrix[j, k] = val > 1e-5 ? val : 0.0
        end
    end
    
    # --- PHYSICAL TIME BOUNDS (Change of Variables!) ---
    safe_u_min = max(0.5, vehicle_params.u_min) 
    safe_u_max = 5.0 
    
    dt_min = ds / safe_u_max
    dt_max = ds / safe_u_min
    
    @variable(model, dt_min <= dt_var[1:N] <= dt_max, start = ds/1.75)
    
    @NLexpression(model, u_expr[k=1:N], ds / dt_var[k])
    @expression(model, T_lap, sum(dt_var[k] for k in 1:N))
    
    @variable(model, 0 <= q[1:N, 1:N+1] <= 1.0, start=0.5)
    @variable(model, node_clarity[1:N] >= 0)
    
    # 1. Calculate the true Time-Integral of clarity for each node
    for j in 1:N
        @NLconstraint(model, 
            node_clarity[j] == sum( q[j, k+1] * dt_var[k] for k in 1:N )
        )
    end
    
    # 2. THE EPIGRAPH VARIABLE
    @variable(model, avg_clarity >= 0)
    
    # 3. THE EPIGRAPH CONSTRAINT
    @NLconstraint(model, 
        avg_clarity * T_lap <= sum(weights[j] * node_clarity[j] for j in 1:N)
    )
    
    # 4. PURELY LINEAR OBJECTIVE
    @objective(model, Max, avg_clarity)
    
    # Energy Constraint
    scale_factor = 1e6
    @NLconstraint(model, 
        sum( (vehicle_params.kh * dt_var[k] + vehicle_params.km * (u_expr[k]^2) * ds) / scale_factor for k in 1:N) <= (P_in_W_avg * T_lap) / scale_factor
    )
    
    # --- EXACT NONLINEAR DYNAMICS (Implicit Euler) ---
    for k in 1:N
        k_next = (k == N) ? 1 : k + 1 
        for j in 1:N
            S_k = S_matrix[j, k]
            S_k_next = S_matrix[j, k_next]
            
            if S_k == 0.0 && S_k_next == 0.0
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] - dt_var[k] * alpha_decay * q[j, k+1]^2
                )
            else
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] + dt_var[k] * (
                        S_k_next * (1 - q[j, k+1])^2 - alpha_decay * q[j, k+1]^2
                    )
                )
            end
        end
    end

    # Slack Periodic Boundaries
    @constraint(model, [j=1:N], q[j, 1] - q[j, N+1] <= 0.01)
    @constraint(model, [j=1:N], q[j, N+1] - q[j, 1] <= 0.01)
    
    optimize!(model)
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]
        println("Success! Speeds optimized.")
        return [ds / v for v in value.(dt_var)] 
    else
        println("Warning: Optimizer failed. Status: ", termination_status(model))
        return fill(1.75, N) 
    end
end

function get_2d_position(s_mod, s_vec, lon_path, lat_path)
    if s_mod >= s_vec[end]
        return lon_path[end], lat_path[end]
    end
    
    idx = findfirst(x -> x >= s_mod, s_vec)
    if idx === nothing || idx == 1
        return lon_path[1], lat_path[1]
    end
    
    s0, s1 = s_vec[idx-1], s_vec[idx]
    frac = (s_mod - s0) / (s1 - s0)
    lon_r = lon_path[idx-1] + frac * (lon_path[idx] - lon_path[idx-1])
    lat_r = lat_path[idx-1] + frac * (lat_path[idx] - lat_path[idx-1])
    return lon_r, lat_r
end

end # module IVESim