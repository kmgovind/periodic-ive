module IVESim

using Ipopt, JuMP, .Clarity

export optimize_lap_speeds_full_spatiotemporal, get_2d_position

function optimize_lap_speeds_full_spatiotemporal(weights, N, ds, P_in_W_avg, alpha_decay, vehicle_params, sigma_val)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # NEW 1: Compute the true spatial integral of the weights
    W_int = sum(weights) * ds 
    
    # Precompute the sensing matrix S[j, k]
    S_matrix = zeros(N, N)
    path_length = N * ds
    for j in 1:N
        for k in 1:N
            dist_direct = abs((j-1)*ds - (k-1)*ds)
            dist_wrap = path_length - dist_direct
            shortest_dist = min(dist_direct, dist_wrap)
            
            # SPARSIFY: Force tiny values to strictly 0.0
            val = exp(-(shortest_dist^2) / (2 * sigma_val^2))
            S_matrix[j, k] = val > 1e-4 ? val : 0.0
        end
    end
    
    safe_u_min = max(0.1, vehicle_params.u_min) 
    
    # NEW 2: Remove the rigid u_max limit to allow sprinting. 
    # (We keep a tiny 1e-4 lower bound just to prevent divide-by-zero errors in Ipopt)
    @variable(model, 1e-4 <= dt[1:N] <= ds/safe_u_min)
    
    @variable(model, 0 <= q[1:N, 1:N+1] <= 1.0)
    
    @variable(model, T_lap >= 0)
    @constraint(model, T_lap == sum(dt[k] for k in 1:N))
    
    # NEW 3: Scale the J_avg limits using the spatial integral
    @variable(model, 0 <= J_avg <= W_int)
    
    # NEW 4: Multiply by `ds` inside the sum to properly integrate across space
    @NLconstraint(model, 
        J_avg * T_lap == sum( weights[j] * 0.5 * (q[j, k] + q[j, k+1]) * dt[k] * ds for j in 1:N, k in 1:N )
    )
    @objective(model, Max, J_avg)
    
    @variable(model, E_avail >= 0)
    @constraint(model, E_avail == P_in_W_avg * T_lap)
    
    # Power Budget (Restored to the proper inequality per your toy problem)
    @NLconstraint(model, 
        sum(vehicle_params.kh * dt[k] + vehicle_params.km * (ds^3 / dt[k]^2) for k in 1:N) <= E_avail
    )
    
    for k in 1:N
        k_next = (k == N) ? 1 : k + 1 
        for j in 1:N
            S_k = S_matrix[j, k]
            S_k_next = S_matrix[j, k_next]
            @NLconstraint(model, 
                q[j, k+1] == q[j, k] + 0.5 * dt[k] * (
                    (S_k * (1 - q[j, k])^2 - alpha_decay * q[j, k]^2) + 
                    (S_k_next * (1 - q[j, k+1])^2 - alpha_decay * q[j, k+1]^2)
                )
            )
        end
    end
    
    @constraint(model, [j=1:N], q[j, 1] == q[j, N+1])
    
    u_nominal = 0.5 * (safe_u_min + vehicle_params.u_max)
    set_start_value.(dt, ds / u_nominal)
    set_start_value(T_lap, N * ds / u_nominal)
    set_start_value.(q, 0.5)
    
    # NEW 5: Warm start scaled by the spatial integral
    set_start_value(J_avg, 0.5 * W_int)
    set_start_value(E_avail, P_in_W_avg * (N * ds / u_nominal))
    
    optimize!(model)
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        return ds ./ value.(dt) 
    else
        println("Warning: Optimizer failed. Defaulting to nominal speed.")
        return fill(vehicle_params.u_max, N) 
    end
end

function get_2d_position(s_mod, s_vec, lon_path, lat_path)
    # Catch floating-point overflow at the end of the lap
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
