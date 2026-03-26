module IVESim

using Ipopt, JuMP

export optimize_lap_speeds_full_spatiotemporal, get_2d_position

function optimize_lap_speeds_full_spatiotemporal(weights, q_initial, N, ds, P_in_W_avg, alpha_decay, vehicle_params, sigma_val)
    model = Model(Ipopt.Optimizer)
    
    # set_silent(model)

    set_attribute(model, "tol", 1e-4)
    set_attribute(model, "acceptable_tol", 1e-3)
    set_attribute(model, "max_iter", 1000) 
    set_attribute(model, "max_cpu_time", 600.0)
    set_attribute(model, "mumps_mem_percent", 10000)
    
    S_matrix = zeros(N, N)
    path_length = N * ds
    for j in 1:N
        for k in 1:N
            dist_direct = abs((j-1)*ds - (k-1)*ds)
            dist_wrap = path_length - dist_direct
            shortest_dist = min(dist_direct, dist_wrap)
            val = 0.005 * exp(-(shortest_dist^2) / (2 * sigma_val^2))
            S_matrix[j, k] = val > 1e-6 ? val : 0.0
        end
    end
    
    # --- PHYSICAL SPEED BOUNDS ---
    # Prevents 1/u derivative explosion while staying strictly positive u > 0
    safe_u_min = max(0.5, vehicle_params.u_min) # Steerage limit
    safe_u_max = 5.0 # Give Ipopt a bounded search space
    
    @variable(model, safe_u_min <= u[1:N] <= safe_u_max, start = 1.75)
    
    @NLexpression(model, dt[k=1:N], ds / u[k])
    @NLexpression(model, T_lap, sum(dt[k] for k in 1:N))
    
    # @variable(model, 0 <= q[1:N, 1:N+1] <= 1.0, start=0.5)
    # Let the nonlinear dynamics handle the upper limit naturally
    @variable(model, q[1:N, 1:N+1] >= 0.0, start=0.5)
    @variable(model, node_clarity[1:N] >= 0)
    
    # --- FIXED TIME INTEGRAL ---
    # Removed '* ds' to properly calculate ∫ q(t) dt
    for j in 1:N
        @NLconstraint(model, 
            node_clarity[j] == sum( 0.5 * (q[j, k] + q[j, k+1]) * dt[k] for k in 1:N )
        )
    end
    
    # Objective exactly matches Eq (2a)
    @NLobjective(model, Max, sum(weights[j] * node_clarity[j] for j in 1:N) / T_lap)
    
    # Energy Constraint exactly matches Eq (2e), scaled by 1e6 for numerical stability
    scale_factor = 1e6
    @NLconstraint(model, 
        sum( (vehicle_params.kh * dt[k] + vehicle_params.km * (u[k]^2) * ds) / scale_factor for k in 1:N) <= (P_in_W_avg * T_lap) / scale_factor
    )
    
    # --- EXACT NONLINEAR DYNAMICS (Eq 2c) ---
    # Restored trapezoidal collocation to match S*(1-q)^2 - alpha*q^2
    for k in 1:N
        k_next = (k == N) ? 1 : k + 1 
        for j in 1:N
            S_k = S_matrix[j, k]
            S_k_next = S_matrix[j, k_next]
            
            if S_k == 0.0 && S_k_next == 0.0
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] - 0.5 * dt[k] * alpha_decay * (q[j, k]^2 + q[j, k+1]^2)
                )
            else
                @NLconstraint(model, 
                    q[j, k+1] == q[j, k] + 0.5 * dt[k] * (
                        (S_k * (1 - q[j, k])^2 - alpha_decay * q[j, k]^2) + 
                        (S_k_next * (1 - q[j, k+1])^2 - alpha_decay * q[j, k+1]^2)
                    )
                )
            end
        end
    end
    
    # REMOVE THIS (Paradigm 2 / MPC):
    # @constraint(model, [j=1:N], q[j, 1] == q_initial[j])

    # ADD THIS (Paradigm 1 / Steady-State Periodic Orbit):
    @constraint(model, [j=1:N], q[j, 1] == q[j, N+1])
    
    optimize!(model)
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]
        println("Success! Speeds optimized.")
        return value.(u) 
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