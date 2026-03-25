# Only include and load packages if they haven't been loaded yet in this session
if !@isdefined(CBOFSData)
    println("Loading modules and packages...")
    include("../src/cbdata.jl")
    include("../src/ASVGeometry.jl")
    include("../src/clarity.jl")
    include("../src/IVESim.jl")
    include("../src/ntfy.jl")

    using .CBOFSData, .ASVGeometry, .Clarity, .IVESim, .Notify
    using Dates, Statistics, Plots, JuMP, Ipopt, ProgressMeter, JLD2
    
    # Optional: If you are actively editing cbdata.jl, ASVGeometry.jl, etc. 
    # you can use Revise to automatically track changes without re-including!
    using Revise
end

# Load salinity data
frames = discover_files("./datafiles")

fr = frames[1]
lon_valid, lat_valid, salt_valid = load_surface_salinity(fr.file)
roi = (lat_valid .>= lat_min) .& (lat_valid .<= lat_max)


# Initialize ASV geometry and simulate path discretization
u_nominal = 1.75  # m/s
vehicle_params = VehicleParams(u_nominal, 24*3600.0)  # 24 hours at 1.75 m/s nominal speed

npts = 200
lon_path, lat_path, s_vec = discretize_polyline(transect_lon_shifted, transect_lat_shifted, npts)
path_length_m = s_vec[end]

# Precompute Local Flat-Earth Cartesian Coordinates for the Path
R_earth = 6371000.0
mean_lat_rad = deg2rad(mean(lat_path))
cos_mean_lat = cos(mean_lat_rad)

# X and Y coordinates in meters
X_path = deg2rad.(lon_path) .* (R_earth * cos_mean_lat)
Y_path = deg2rad.(lat_path) .* R_earth


println("running optimizer debugging....")
# Simulation parameters
dt_step_sec = 60.0              # 1 minute per step
sigma_val = 300.0              # Sensing footprint [m]
alpha_val = 0.001                # Clarity decay

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
# current_weights = fill(1.0, npts)

gaussian_weights(x, mu, sigma) = exp(-((x - mu) / sigma)^2)
weights_x = range(0, 1, length=npts)
current_weights = gaussian_weights.(weights_x, 0.5, 0.15)

q_initial = fill(0.0, npts)  # Initial clarity guess for optimization


# Run optimizer
opt_time_start = time()
u_opt_segments = optimize_lap_speeds_full_spatiotemporal(current_weights, q_initial, N_segments, ds, P_in_W, alpha_val, vehicle_params, sigma_val)
opt_time_elapsed = time() - opt_time_start


send_ntfy("Debugging completed in $(round(opt_time_elapsed, digits=2)) seconds", "Sim Finished", "default");

plot(u_opt_segments, title="Optimized Speeds Along Path", xlabel="Segment Index", ylabel="Speed (m/s)")
plot!(current_weights[1:end-1], title="Optimized Speeds Along Path", xlabel="Segment Index", ylabel="Speed (m/s)", label="Weights", color=:red)
