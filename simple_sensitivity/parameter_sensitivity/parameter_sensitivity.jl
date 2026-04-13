using JuMP, Ipopt, Plots, QuadGK, DataFrames, Statistics, JLD2, Dates

# ==========================================
# 1. Plot Settings & Global Constants
# ==========================================
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2,
    guidefontsize=14,
    tickfontsize=14,
    legendfontsize=10,
    titlefontsize=16
)

const L = 100.0             # Total Path Length (m)
const P_budget = 750.0      # Average Power In (W)
const u_min = 0.01           # Min feasible speed (m/s)
const u_max_guess = 2.5     # Used for warm starts
const kh = 10.0             # Hotel Load (W)
const km = 83.0             # Drag Coefficient
const S_scale = 10.0        # Sensing magnitude

# ==========================================
# 2. Environment & Physics Helpers
# ==========================================
function sensing_rate_continuous(s)
    peak1 = 0.5 * exp(-((s - 30)^2) / (2 * 5^2))
    peak2 = 0.8 * exp(-((s - 80)^2) / (2 * 2^2))
    raw = 0.05 + S_scale * (peak1 + peak2)
    raw_min = 0.05
    raw_max = 0.05 + S_scale * 0.8  
    # Clamp the bottom to 0.01 to prevent singular Jacobians
    return clamp((raw - raw_min) / (raw_max - raw_min), 0.01, 1.0)
end

function get_segment_S(s_start, s_end)
    val, err = quadgk(sensing_rate_continuous, s_start, s_end)
    return val / (s_end - s_start)
end

function cyclic_distance(s1, s2, L)
    dist = abs(s1 - s2)
    return min(dist, L - dist)
end

function get_footprint_matrix(N, L; sigma=10.0)
    ds = L / N
    footprint = zeros(N, N)
    for i in 1:N
        s_robot = (i - 0.5) * ds
        for j in 1:N
            s_target = (j - 0.5) * ds
            dist = cyclic_distance(s_robot, s_target, L)
            footprint[i, j] = exp(-(dist^2) / (2 * sigma^2))
        end
    end
    return footprint
end

# ==========================================
# 3. Parameterized Optimizer
# ==========================================
# Notice alpha and sigma are now keyword arguments!
function solve_optimal_trajectory(N; P_in=750.0, alpha=0.1, sigma=10.0)
    ds = L / N
    gamma_vals = [get_segment_S((i-1)*ds, i*ds) for i in 1:N]
    Gamma_total = sum(gamma_vals) * ds 
    
    S_base = 1.0
    footprint_matrix = get_footprint_matrix(N, L, sigma=sigma) 
    
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # Floor dt at ds/10.0 to prevent 1/dt^3 singularities
    @variable(model, ds/10.0 <= dt[1:N] <= ds/u_min)  
    @variable(model, 0 <= q[1:N+1, 1:N] <= 1.0)
    
    @variable(model, T_lap >= 0)
    @constraint(model, T_lap == sum(dt[i] for i in 1:N))
    
    @variable(model, 0 <= J_avg <= Gamma_total)
    
    # Fast O(N) Spatial Integration
    @variable(model, spatial_integral[1:N+1] >= 0)
    @constraint(model, [i=1:N+1], 
        spatial_integral[i] == sum(gamma_vals[j] * q[i, j] * ds for j in 1:N)
    )
    @NLconstraint(model, J_avg * T_lap == sum( 0.5 * (spatial_integral[i] + spatial_integral[i+1]) * dt[i] for i in 1:N ))
    @objective(model, Max, J_avg)
    
    # Energy
    @variable(model, E_avail >= 0)
    @constraint(model, E_avail == P_in * T_lap)
    @NLconstraint(model, sum(kh * dt[i] + km * (ds / dt[i])^3 * dt[i] for i in 1:N) <= E_avail)
    
    # Clarity Dynamics parameterized by local `alpha` and `sigma` footprint
    for i in 1:N          
        for j in 1:N      
            S_val = S_base * footprint_matrix[i, j]
            @NLconstraint(model, 
                q[i+1, j] == q[i, j] + 0.5 * dt[i] * (
                    (S_val * (1 - q[i, j])^2 - alpha * q[i, j]^2) + 
                    (S_val * (1 - q[i+1, j])^2 - alpha * q[i+1, j]^2)
                )
            )
        end
    end
    
    @constraint(model, [j=1:N], q[1, j] == q[N+1, j])
    
    # Warm Start
    u_nominal = 0.5 * (u_min + u_max_guess)
    set_start_value.(dt, ds / u_nominal)
    set_start_value.(q, 0.5)                 
    set_start_value(T_lap, N * ds / u_nominal)
    set_start_value(J_avg, 0.5 * Gamma_total)  
    set_start_value.(spatial_integral, 0.5 * Gamma_total)
    set_start_value(E_avail, P_in * (N * ds / u_nominal))

    # Solve
    t_start = time()
    optimize!(model)
    solve_time = time() - t_start
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        return objective_value(model), solve_time, value.(dt)
    else
        println("Warning: Optimizer failed for alpha=$alpha, sigma=$sigma")
        return NaN, solve_time, fill(NaN, N)
    end
end

# ==========================================
# 4. Sensitivity Sweep Runner
# ==========================================
N_sweep = 100 # Keep N relatively low for the sweep to run quickly

alpha_vals = [0.01, 0.05, 0.1, 0.2, 0.5]
sigma_vals = [2.0, 5.0, 10.0, 20.0, 50.0]

num_a = length(alpha_vals)
num_s = length(sigma_vals)

obj_matrix = zeros(num_a, num_s)
time_matrix = zeros(num_a, num_s)
trajectories_dict = Dict()

println("Starting 2D Sensitivity Sweep (N=$N_sweep)...")
for (i, a) in enumerate(alpha_vals)
    for (j, s) in enumerate(sigma_vals)
        print("Solving α = $a, σ = $s ... ")
        obj, t, dts = solve_optimal_trajectory(N_sweep, P_in=P_budget, alpha=a, sigma=s)
        
        obj_matrix[i, j] = obj
        time_matrix[i, j] = t
        trajectories_dict[(a, s)] = (L/N_sweep) ./ dts # Store speed profile
        
        println("Obj = $(round(obj, digits=2)), Time = $(round(t, digits=2))s")
    end
end

# Save the data to a jld2 file for later analysis
datetime_str = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
filename = datetime_str * "_sensitivity_results.jld2"
mkdir("$datetime_str")
filepath = joinpath("$datetime_str", filename)
@save filepath alpha_vals sigma_vals obj_matrix time_matrix trajectories_dict

# ==========================================
# 5. Plotting the Results
# ==========================================
# A. Heatmap: Objective Value
p_obj = heatmap(string.(sigma_vals), string.(alpha_vals), obj_matrix, 
    xlabel="Sensing Radius (σ) [m]", ylabel="Decay Rate (α)", 
    title="Objective Value Sensitivity", color=:viridis,
    right_margin=5Plots.mm)
savefig(p_obj, "$datetime_str/sensitivity_objective.png")

# B. Heatmap: Compute Time
p_time = heatmap(string.(sigma_vals), string.(alpha_vals), time_matrix, 
    xlabel="Sensing Radius (σ) [m]", ylabel="Decay Rate (α)", 
    title="Compute Time Sensitivity (s)", color=:inferno,
    right_margin=5Plots.mm)
savefig(p_time, "$datetime_str/sensitivity_time.png")

# C. Speed Profiles: Varying Alpha (Fixed Sigma = 10.0)
p_alpha_sweep = plot(title="Speed Profiles vs Decay Rate (σ = 10.0)", 
    xlabel="Distance (m)", ylabel="Speed (m/s)", legend=:topleft)
s_plot = range(0, L, length=N_sweep)
for a in alpha_vals
    speeds = trajectories_dict[(a, 10.0)]
    plot!(p_alpha_sweep, s_plot, speeds, linetype=:step, lw=2, label="α = $a")
end
savefig(p_alpha_sweep, "$datetime_str/speed_vs_alpha.png")

# D. Speed Profiles: Varying Sigma (Fixed Alpha = 0.1)
p_sigma_sweep = plot(title="Speed Profiles vs Sensing Radius (α = 0.1)", 
    xlabel="Distance (m)", ylabel="Speed (m/s)", legend=:topleft)
for s in sigma_vals
    speeds = trajectories_dict[(0.1, s)]
    plot!(p_sigma_sweep, s_plot, speeds, linetype=:step, lw=2, label="σ = $s")
end
savefig(p_sigma_sweep, "$datetime_str/speed_vs_sigma.png")

println("Sweep complete! Heatmaps and speed profiles saved.")


# ==========================================
# 6. Constant Speed Baseline Comparison
# ==========================================
println("\nStarting Constant Speed Baseline Evaluation...")

# 1. Calculate the energy-limited constant speed
const u_const = ((P_budget - kh) / km)^(1/3)
println("Calculated Constant Speed Baseline: $(round(u_const, digits=3)) m/s")

# 2. Forward simulation function for constant speed
function evaluate_constant_speed(alpha, sigma; N_eval=100)
    ds = L / N_eval
    gamma_vals = [get_segment_S((i-1)*ds, i*ds) for i in 1:N_eval]
    footprint_matrix = get_footprint_matrix(N_eval, L, sigma=sigma)
    
    q_curr = fill(0.5, N_eval)
    q_next = zeros(N_eval)
    
    dt = ds / u_const
    T_lap = L / u_const
    
    num_laps = 15 # Run enough laps to reach steady state
    final_score = 0.0
    
    for lap in 1:num_laps
        total_clarity = 0.0
        for i in 1:N_eval
            for j in 1:N_eval
                S_sensor = 1.0 * footprint_matrix[i, j]
                f(q_val) = S_sensor * (1 - q_val)^2 - alpha * q_val^2
                
                # RK4 Step
                k1 = dt * f(q_curr[j])
                k2 = dt * f(q_curr[j] + k1/2)
                k3 = dt * f(q_curr[j] + k2/2)
                k4 = dt * f(q_curr[j] + k3)
                q_next[j] = q_curr[j] + (k1 + 2k2 + 2k3 + k4)/6
                
                if lap == num_laps
                    total_clarity += 0.5 * gamma_vals[j] * (q_curr[j] + q_next[j]) * dt * ds
                end
            end
            q_curr .= q_next
        end
        if lap == num_laps
            final_score = total_clarity / T_lap
        end
    end
    
    return final_score
end

# 3. Evaluate the grid
const_obj_matrix = zeros(num_a, num_s)
improvement_matrix = zeros(num_a, num_s)

for (i, a) in enumerate(alpha_vals)
    for (j, s) in enumerate(sigma_vals)
        base_obj = evaluate_constant_speed(a, s, N_eval=N_sweep)
        const_obj_matrix[i, j] = base_obj
        
        # Calculate percentage improvement
        opt_obj = obj_matrix[i, j]
        improvement_matrix[i, j] = ((opt_obj - base_obj) / base_obj) * 100.0
    end
end

# Save the new baseline data to the JLD2 file
@save filepath alpha_vals sigma_vals obj_matrix time_matrix trajectories_dict const_obj_matrix improvement_matrix

# 4. Plot the Percentage Improvement Heatmap
p_imp = heatmap(string.(sigma_vals), string.(alpha_vals), improvement_matrix, 
    xlabel="Sensing Radius (σ) [m]", ylabel="Decay Rate (α)", 
    title="Optimal Controller Improvement (%)", color=:cividis,
    right_margin=5Plots.mm)

# Add text annotations on top of the heatmap blocks so the exact % is readable
for i in 1:num_a
    for j in 1:num_s
        val = round(improvement_matrix[i, j], digits=1)
        annotate!(p_imp, j, i, text("$(val)%", 10, :white, :center))
    end
end

savefig(p_imp, "$datetime_str/sensitivity_improvement.png")
println("Baseline evaluation complete! Improvement heatmap saved.")