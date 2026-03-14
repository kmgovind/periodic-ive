module CBOFSData

using NCDatasets
using Dates
using Printf

export CBOFSFrame, discover_files, load_surface_salinity

# A simple struct to keep your metadata organized
struct CBOFSFrame
    file::String
    dt::DateTime
end

"""
    discover_files(directory=".")
Scans a directory for CBOFS NetCDF files and returns a sorted list of CBOFSFrame objects.
"""
function discover_files(path=".")
    pat = r"^chesapeake_salinity_(\d{8})_t(\d{2})z_n(\d{3})\.nc$"
    all_files = readdir(path)
    
    frames = CBOFSFrame[]
    
    for f in all_files
        m = match(pat, f)
        if m !== nothing
            ymd = m.captures[1]
            cyc = parse(Int, m.captures[2])
            fh  = parse(Int, m.captures[3])
            
            # Construct exact UTC timestamp
            dt = DateTime(ymd, dateformat"yyyymmdd") + Hour(cyc) + Hour(fh)
            push!(frames, CBOFSFrame(joinpath(path, f), dt))
        end
    end
    
    # Sort by datetime so your simulation proceeds chronologically
    sort!(frames, by = x -> x.dt)
    return frames
end

"""
    load_surface_salinity(file_path)
Extracts 1D vectors of Lon, Lat, and Salinity for the surface layer.
"""
function load_surface_salinity(file_path::String)
    if !isfile(file_path)
        error("File not found: $file_path")
    end

    ds = NCDataset(file_path)
    try
        # Using [:,:] ensures we get the full 2D grid
        lon = ds["lon_rho"][:,:]
        lat = ds["lat_rho"][:,:]
        # Surface layer is the last index of the 3rd dimension (s_rho)
        salt = ds["salt"][:, :, end, 1]
        
        # Create mask: Filter missing, NaN, and land-fill values
        # ROMS typically uses 1e37 for land; salinity > 100 is impossible
        mask = .!ismissing.(salt) .&& .!isnan.(salt) .&& (salt .< 100.0)
        
        # Convert to 1D vectors for plotting/interpolation
        lon_v  = vec(lon)[vec(mask)]
        lat_v  = vec(lat)[vec(mask)]
        salt_v = Float64.(vec(salt)[vec(mask)])
        
        return lon_v, lat_v, salt_v
    finally
        close(ds) # Always close in a 'finally' block to avoid locked files
    end
end

end # module