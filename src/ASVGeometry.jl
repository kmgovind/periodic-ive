module ASVGeometry

using LinearAlgebra, QuadGK

export discretize_polyline, sample_salinity
export VehicleParams, SimParams

Base.@kwdef struct VehicleParams
    u_min::Float64 = 0.75
    u_max::Float64 = 2.5
    kh::Float64 = 10.0
    km::Float64 = 83.0
    E_budget::Float64 = 0.0
end

function VehicleParams(u_nominal::Float64, duration_sec::Float64; u_min::Float64=0.00, u_max::Float64=2.5, kh::Float64=10.0, km::Float64=83.0)
    # Energy = Power * Time; Power = kh + km * u^3
    # For constant speed u_nominal over duration_sec, compute energy budget in Wh
    E_budget = (kh + km * u_nominal^3) * duration_sec / 3600.0  # Convert Ws to Wh
    return VehicleParams(u_min=u_min, u_max=u_max, kh=kh, km=km, E_budget=E_budget)
end
    

Base.@kwdef struct SimParams
    alpha::Float64 = 0.05
    S_scale::Float64 = 1.0
    sigma::Float64 = 1000.0 # Standard deviation for sensing dropoff in meters
end

function discretize_polyline(lon::Vector{<:Real}, lat::Vector{<:Real}, npts::Int)
    n = length(lon)
    @assert n == length(lat) "lon/lat must have same length"
    @assert n >= 2 "need at least 2 points"

    # Flat earth approximation for fast metric distance
    φ0 = deg2rad(sum(lat)/length(lat))
    m_per_deg_lat = 111_320.0
    m_per_deg_lon = 111_320.0 * cos(φ0)

    seglen = zeros(Float64, n - 1)
    for i in 1:n-1
        dx = (lon[i+1] - lon[i]) * m_per_deg_lon
        dy = (lat[i+1] - lat[i]) * m_per_deg_lat
        seglen[i] = hypot(dx, dy)
    end

    s = vcat(0.0, cumsum(seglen))
    st = range(0.0, s[end], length=npts)

    lon_out = zeros(Float64, npts)
    lat_out = zeros(Float64, npts)

    j = 1
    for k in eachindex(st)
        while j < length(s)-1 && st[k] > s[j+1]
            j += 1
        end
        ds = s[j+1] - s[j]
        ξ = ds > 0 ? (st[k] - s[j]) / ds : 0.0
        lon_out[k] = (1-ξ) * lon[j] + ξ * lon[j+1]
        lat_out[k] = (1-ξ) * lat[j] + ξ * lat[j+1]
    end

    # Return coordinates and the cumulative distance vector
    return lon_out, lat_out, collect(st)
end

function sample_salinity(lon_r, lat_r, lon_v, lat_v, salt_v)
    φ0 = deg2rad(lat_r)
    m_per_deg_lat = 111_320.0
    m_per_deg_lon = 111_320.0 * cos(φ0)
    
    # Squared distance is faster than hypot/sqrt for nearest-neighbor checks
    d2 = ((lon_v .- lon_r) .* m_per_deg_lon).^2 .+ ((lat_v .- lat_r) .* m_per_deg_lat).^2
    return salt_v[argmin(d2)]
end

end # module ASVGeometry