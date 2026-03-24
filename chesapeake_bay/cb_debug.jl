include("../src/cbdata.jl")
include("../src/ASVGeometry.jl")
include("../src/clarity.jl")

using .CBOFSData, .ASVGeometry, .Clarity
using Dates, Statistics, Plots, JuMP, Ipopt, Revise, ProgressMeter


############# Optimizer ###########################################################################################################################
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
            S_matrix[j, k] = exp(-(shortest_dist^2) / (2 * sigma_val^2))
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
########################################################################################################################################################


################### Transect Definition #######################################################################################
# Start at the southeastern-most corner; southern return is offset by the same spacing (0.05°) as leg spacing
transect_lon = [
    -75.95, -75.95, -75.95,  # Start SE corner -> bottom of Leg 1 -> top of Leg 1
    -76.00, -76.00,          # Across top -> down Leg 2
    -76.05, -76.05,          # Across bottom -> up Leg 3
    -76.10, -76.10,          # Across top -> down Leg 4
    -76.15, -76.15,          # Across bottom -> up Leg 5
    -76.20, -76.20, -76.20,  # Across top -> down Leg 6 -> drop to southern return
    -75.95                   # Return east to close loop at start
]

t_lat_max = 37.05;
t_lat_min = 36.90;
t_lat_return = 36.875;  # Southern return line latitude

transect_lat = [
    t_lat_return, t_lat_min, t_lat_max,  # Start SE -> Leg 1 bottom -> Leg 1 top
    t_lat_max, t_lat_min,         # Leg 2 top -> Leg 2 bottom
    t_lat_min, t_lat_max,         # Leg 3 bottom -> Leg 3 top
    t_lat_max, t_lat_min,         # Leg 4 top -> Leg 4 bottom
    t_lat_min, t_lat_max,         # Leg 5 bottom -> Leg 5 top
    t_lat_max, t_lat_min, t_lat_return,  # Leg 6 top -> Leg 6 bottom -> southern return line
    t_lat_return                 # Close loop at SE start
]

# Shift path slightly west and north to avoid land overlap
lon_shift = -0.065   # west
lat_shift =  0.1  # north
transect_lon_shifted = transect_lon .+ lon_shift
transect_lat_shifted = transect_lat .+ lat_shift

lat_min, lat_max = 36.75, 37.5

### Load salinity data
frames = discover_files("./datafiles")

fr = frames[1]
lon_valid, lat_valid, salt_valid = load_surface_salinity(fr.file)
roi = (lat_valid .>= lat_min) .& (lat_valid .<= lat_max)


###############################  Simulation   ###########################################

# Initialize ASV geometry and simulate path discretization
u_nominal = 1.75  # m/s
vehicle_params = VehicleParams(u_nominal, 24*3600.0)  # 24 hours at 1.75 m/s nominal speed

npts = 200
lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
path_length_m = s_vec[end]


######################### Optimizer Debugging ################
println("running optimizer debugging....")
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
current_weights = fill(1.0, npts)


# Run optimizer
u_opt_segments_full = optimize_lap_speeds_full_spatiotemporal(current_weights, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val)
plot(u_opt_segments_full, ylims=(0, vehicle_params.u_max * 1.2), xlabel="Segment Index", ylabel="Optimized Speed (m/s)", title="Optimized Segment Speeds (Full Spatiotemporal Model)", legend=false)
