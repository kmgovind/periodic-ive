module Clarity
export clarity_dynamics, sensing_function, update_clarity_rk4, calculate_target_weights

# Dynamics strictly match Equation 2 from the main body text:
function clarity_dynamics(qj, Sj; alpha = 0.05)
    return Sj * (1.0 - qj)^2 - alpha * qj^2
end

# Now accepts true physical distance rather than a generic vector norm
# Matches Equation 3:
function sensing_function(distance_meters; S_0 = 1.0, sigma = 1000.0)
    return S_0 * exp(-(distance_meters^2) / (2 * sigma^2))
end

# Upgraded ODE Integrator (Runge-Kutta 4th Order)
# Provides smooth, accurate integration without needing clamp() hacks
function update_clarity_rk4(qj, Sj, dt; alpha = 0.05)
    k1 = clarity_dynamics(qj, Sj; alpha=alpha)
    k2 = clarity_dynamics(qj + 0.5 * dt * k1, Sj; alpha=alpha)
    k3 = clarity_dynamics(qj + 0.5 * dt * k2, Sj; alpha=alpha)
    k4 = clarity_dynamics(qj + dt * k3, Sj; alpha=alpha)
    
    q_new = qj + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    # Optional failsafe for floating point noise, but RK4 shouldn't overshoot
    return clamp(q_new, 0.0, 1.0) 
end

function calculate_target_weights(measurements_buffer, pos_buffer, s_vec, npts; rated_salinity=30.0, sigma_weight=3.0)

    # Gaussian weight
    weight_exp(s, s_ref=rated_salinity, σ=sigma_weight) = exp(-((s - s_ref)^2) / (2σ^2))

    # Cutoff weight
    weight_ge30(s; threshold=rated_salinity, σ=sigma_weight) = s >= threshold ? 1.0 : exp(-((s - threshold)^2) / (2σ^2))

    new_weights = ones(Float64, npts)
    
    # If no measurements yet (e.g., Lap 1), return uniform weights
    if isempty(measurements_buffer)
        return new_weights
    end
    
    # For each path point, find the closest measurement from the previous lap
    for i in 1:npts
        closest_idx = argmin(abs.(pos_buffer .- s_vec[i]))
        s_meas = measurements_buffer[closest_idx]
        
        # Apply custom weighting function
        new_weights[i] = weight_ge30(s_meas)
    end
    
    # --- CRITICAL SOLVER STEP: Normalization ---
    avg_weight = sum(new_weights) / npts
    if avg_weight > 1e-6 
        new_weights ./= avg_weight
    end
    
    return new_weights
end

end # module Clarity