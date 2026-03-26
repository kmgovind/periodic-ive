using JuMP, Ipopt, Plots, QuadGK, DataFrames, Statistics

# Plot Settings
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2,
    guidefontsize=14,
    tickfontsize=14,
    legendfontsize=10,
    titlefontsize=18
)

# ==========================================
# 1. Physics & Environment
# ==========================================
const L = 100.0             # Total Path Length (m)
const P_budget = 750.0     # Average Power In (W)
const u_min = 0.0           # Min feasible speed (m/s)
const u_max = 2.5           # Max speed (m/s)
const kh = 10.0             # Hotel Load (W)
const km = 83.0              # Drag Coefficient
const alpha = 0.1          # Decay rate
const S_scale = 10.0         # Sensing magnitude

# Environment Map (Gaussian Hotspots)
function sensing_rate_continuous(s)
    peak1 = 0.5 * exp(-((s - 30)^2) / (2 * 5^2))
    peak2 = 0.8 * exp(-((s - 80)^2) / (2 * 2^2)) # Sharper peak

    raw = 0.05 + S_scale * (peak1 + peak2)

    # return raw

    # Normalize to [0, 1]
    raw_min = 0.05
    raw_max = 0.05 + S_scale * 0.8  # dominant peak scale
    return clamp((raw - raw_min) / (raw_max - raw_min), 0.01, 1.0)
end

# function sensing_rate_continuous(s)
#     return 1.0  # Uniform sensing rate for simplicity
# end

# Distance on a periodic 1D track of length L
function cyclic_distance(s1, s2, L)
    dist = abs(s1 - s2)
    return min(dist, L - dist)
end

# Gaussian footprint: S_max is peak rate, sigma is physical spread in meters
function get_S_matrix(N, L, S_max, sigma)
    ds = L / N
    S_mat = zeros(N, N)
    for i in 1:N # Robot location index
        s_robot = (i - 0.5) * ds
        for j in 1:N # Target segment index
            s_target = (j - 0.5) * ds
            dist = cyclic_distance(s_robot, s_target, L)
            S_mat[i, j] = S_max * exp(-(dist^2) / (2 * sigma^2))
        end
    end
    return S_mat
end

# Anti-aliasing: Integrate S over the segment
function get_segment_S(s_start, s_end)
    val, err = quadgk(sensing_rate_continuous, s_start, s_end)
    return val / (s_end - s_start)
end


# Computes the shortest distance on a periodic 1D track of length L
function cyclic_distance(s1, s2, L)
    dist = abs(s1 - s2)
    return min(dist, L - dist)
end

# Generates an NxN matrix of sensor attenuation based on physical distance
function get_footprint_matrix(N, L; sigma=10.0)
    ds = L / N
    footprint = zeros(N, N)
    for i in 1:N # Robot location index
        s_robot = (i - 0.5) * ds
        for j in 1:N # Target segment index
            s_target = (j - 0.5) * ds
            dist = cyclic_distance(s_robot, s_target, L)
            # Gaussian dropoff (1.0 at center, decaying with distance)
            footprint[i, j] = exp(-(dist^2) / (2 * sigma^2))
        end
    end
    return footprint
end

# ==========================================
# Optimization Routine (Time-Averaged Total Clarity)
# ==========================================
function solve_optimal_trajectory(N; P_in=750.0)
    ds = L / N
    S_vals = [get_segment_S((i-1)*ds, i*ds) for i in 1:N]
    
    # NEW: Precompute the sensor footprint matrix
    footprint_matrix = get_footprint_matrix(N, L, sigma=10.0) 
    
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    
    # --- Variables ---
    # @variable(model, ds/u_max <= dt[1:N] <= ds/u_min)
    @variable(model, dt[1:N] <= ds / u_min)  # Only upper bound, since dt can be large if u is small, no more u_max constraint
    @variable(model, 0 <= q[1:N+1, 1:N] <= 1.0)
    
    @variable(model, T_lap >= 0)
    @constraint(model, T_lap == sum(dt[i] for i in 1:N))
    
    # NEW: Bounds for J_avg are now [0, L] because it is the total clarity over length L
    @variable(model, 0 <= J_avg <= L)
    
    # NEW: Multiply by `ds` inside the sum for true spatial integration
    @NLconstraint(model, J_avg * T_lap == sum( 0.5 * (q[i, j] + q[i+1, j]) * dt[i] * ds for i in 1:N, j in 1:N ))
    
    # Linear Objective
    @objective(model, Max, J_avg)
    
    # --- Constraints ---
    @variable(model, E_avail >= 0)
    @constraint(model, E_avail == P_in * T_lap)
    
    @NLconstraint(model, sum(kh * dt[i] + km * (ds / dt[i])^3 * dt[i] for i in 1:N) <= E_avail)
    
    # Clarity Dynamics
    for i in 1:N          # Time step index (Robot is physically at segment i)
        for j in 1:N      # Spatial segment index being updated
            
            # NEW: S is the environment's base rate at j, scaled by the robot's footprint from i to j
            S = S_vals[j] * footprint_matrix[i, j]
            
            @NLconstraint(model, 
                q[i+1, j] == q[i, j] + 0.5 * dt[i] * (
                    (S * (1 - q[i, j])^2 - alpha * q[i, j]^2) + 
                    (S * (1 - q[i+1, j])^2 - alpha * q[i+1, j]^2)
                )
            )
        end
    end
    
    # Cyclic boundary condition for all spatial segments
    @constraint(model, [j=1:N], q[1, j] == q[N+1, j])
    
    # --- Warm Start ---
    u_nominal = 0.5 * (u_min + u_max)
    set_start_value.(dt, ds / u_nominal)
    set_start_value.(q, 0.5)                 
    set_start_value(T_lap, N * ds / u_nominal)
    set_start_value(J_avg, 0.5 * L)          # NEW: Warm start scaled by L
    set_start_value(E_avail, P_in * (N * ds / u_nominal))

    # --- Solve ---
    t_start = time()
    optimize!(model)
    solve_time = time() - t_start
    
    if termination_status(model) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
        return objective_value(model), solve_time, value.(dt), value.(q)
    else
        println("Warning: Optimizer failed.")
        return NaN, NaN, [], []
    end
end


# N_values = [10, 25, 50, 100, 200, 400, 600, 1000, 2500, 5000, 7500, 10000]
# N_values = [10, 25, 50, 100, 200, 400, 600, 1000, 2500]
N_values = [10, 50, 100, 200, 400, 800]
results_obj = []
results_time = []
trajectories = []
sim_times = []
results_q = []

println("Running Grid Independence Sweep...")
for N in N_values
    print("  N = $N ... ")
    obj, t, dts, qs = solve_optimal_trajectory(N, P_in=P_budget)
    # obj, t, dts, qs = old_solve_optimal_trajectory(N)
    push!(results_obj, obj)
    push!(results_time, t)
    push!(sim_times, dts)
    push!(results_q, qs)
    
    # Store speed profile for plotting
    speeds = (L/N) ./ dts
    push!(trajectories, speeds)
    println("Done. J = $(round(obj, digits=2)), Time = $(round(t, digits=3))s")
end

# Create results table
results_table = DataFrame(
    N = N_values,
    ObjectiveValue = round.(results_obj, digits=2),
    ComputeTime = round.(results_time, digits=3)
)

println("\n" * "="^60)
println("Grid Independence Sweep Results")
println("="^60)
println(results_table)
println("="^60)

# Export to LaTeX
open("results_table.tex", "w") do f
    write(f, "\\begin{table}\n")
    write(f, "    \\centering\n")
    write(f, "    \\begin{tabular}{ccc}\n")
    write(f, "        \\hline\n")
    write(f, "        N & Objective Function Value & Compute Time (s) \\\\\n")
    write(f, "        \\hline\n")
    for i in eachindex(N_values)
        write(f, "        $(N_values[i]) & $(round(results_obj[i], digits=2)) & $(round(results_time[i], digits=3))\\\\\n")
    end
    write(f, "        \\hline\n")
    write(f, "    \\end{tabular}\n")
    write(f, "    \\caption{Objective function values and computation times for varying grid resolutions.}\n")
    write(f, "    \\label{tab:grid_resolution}\n")
    write(f, "\\end{table}\n")
end

println("LaTeX table exported to results_table.tex")


########## Compute Time Plots ###############
p2 = plot(N_values, results_time, marker=:square, color=:red,
    xaxis=:log, yaxis=:log, title="Compute Time",
    ylabel="Time (s)", xlabel="N", lw=2, label="Simulation Time")
idx = findfirst(!isnan, results_time)
if idx === nothing
    k1 = 1.0
    k2 = 1.0
    k3 = 1.0
else
    k1 = results_time[idx] / N_values[idx]
    k2 = results_time[idx] / (N_values[idx]^2)
    k3 = results_time[idx] / (N_values[idx] * log(N_values[idx]))
    k4 = results_time[idx] / (N_values[idx] * log(N_values[idx])^2)
end

plot!(p2, N_values, k1 .* N_values, lw=2, ls=:dash, color=:green, label="O(N)")
plot!(p2, N_values, k2 .* (N_values .^ 2), lw=2, ls=:dot, color=:blue, label="O(N^2)")
plot!(p2, N_values, k3 .* N_values .* log.(N_values), lw=2, ls=:dashdot, color=:purple, label="O(N log N)")
plot!(p2, N_values, k4 .* N_values .* log.(N_values).^2, lw=2, ls=:dot, color=:orange, label="O(N log^2 N)")

plot!(legend=:topleft)

savefig(p2, "compute_time.png")
##############################################





############# Resolution Comparison Plots ###############
# Speed Profile Comparison
idx_coarse = findfirst(==(50), N_values)
idx_mid = findfirst(==(200), N_values)
idx_fine = findfirst(==(800), N_values)

p3 = plot(
    title = "Resolution Comparison",
    xlabel = "Distance (m)",
    ylabel = "Speed (m/s)",
    ylims = (0, u_max + 1.0),
    legend=:topleft
)

# Left axis: speed profiles
s_coarse = range(0, L, length=length(trajectories[idx_coarse]))
plot!(p3, s_coarse, trajectories[idx_coarse], linetype=:step, label="N=$(N_values[idx_coarse])", lw=2)

s_mid = range(0, L, length=length(trajectories[idx_mid]))
plot!(p3, s_mid, trajectories[idx_mid], linetype=:step, label="N=$(N_values[idx_mid])", lw=2)

s_fine = range(0, L, length=length(trajectories[idx_fine]))
plot!(p3, s_fine, trajectories[idx_fine], linetype=:step, label="N=$(N_values[idx_fine])", lw=2, color=:black)

# hline!(p3, [u_max], lw=2, ls=:dash, color=:red, label="Max Speed = $(u_max) m/s")
# ylims!(p3, 0, u_max + 0.5)

# Right axis: info map
p3r = twinx(p3)
s_map = 0:0.5:L
plot!(
    p3r,
    s_map,
    sensing_rate_continuous.(s_map),
    color = :blue,
    fill = (0, 0.2, :blue),
    alpha = 0.3,
    ylabel = "Info Map Weight",
    label = "Info Map",
    lw = 2,
    ylims=(0.0, 1.0),
    legend=:topright
)

plot!(dpi=1200)
savefig(p3, "grid_independence_sweep.pdf")
#########################################################################