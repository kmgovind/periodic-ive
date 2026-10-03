#!/usr/bin/env julia
# Download the exact CBOFS surface-salinity input window used by the ACC paper.
#
# Usage (from code/chesapeake_bay):
#   julia --project=.. download_cbofs_data.jl
#
# Optional environment overrides:
#   START_DATE=2026-02-16 END_DATE=2026-02-20 DATA_DIR=datafiles \
#     julia --project=.. download_cbofs_data.jl

using Dates
using Downloads

const DEFAULT_START_DATE = Date(2026, 2, 16)
const DEFAULT_END_DATE = Date(2026, 2, 20)
const CYCLES = ("00", "06", "12", "18")

parse_date_env(name, default) = Date(get(ENV, name, string(default)))

function available_forecast_hours(date, cycle)
    yyyy = Dates.format(date, "yyyy")
    mm = Dates.format(date, "mm")
    dd = Dates.format(date, "dd")
    yyyymmdd = Dates.format(date, "yyyymmdd")
    catalog_url = "https://opendap.co-ops.nos.noaa.gov/thredds/catalog/NOAA/CBOFS/MODELS/$yyyy/$mm/$dd/catalog.html"
    html_buffer = IOBuffer()
    Downloads.download(catalog_url, html_buffer)
    html = String(take!(html_buffer))
    pattern = Regex("cbofs\\.t$(cycle)z\\.$yyyymmdd\\.fields\\.n(\\d{3})\\.nc")
    return sort!(unique(match.captures[1] for match in eachmatch(pattern, html)))
end

function download_window(start_date, end_date, output_dir)
    start_date <= end_date || error("START_DATE must not be after END_DATE")
    mkpath(output_dir)
    downloaded = 0
    skipped = 0
    failed = 0

    for date in start_date:Day(1):end_date
        yyyy = Dates.format(date, "yyyy")
        mm = Dates.format(date, "mm")
        dd = Dates.format(date, "dd")
        yyyymmdd = Dates.format(date, "yyyymmdd")
        for cycle in CYCLES
            hours = try
                available_forecast_hours(date, cycle)
            catch err
                @warn "Could not read the NOAA CBOFS catalog" date cycle exception=(err, catch_backtrace())
                failed += 1
                continue
            end
            for hour in hours
                filename = "chesapeake_salinity_$(yyyymmdd)_t$(cycle)z_n$(hour).nc"
                destination = joinpath(output_dir, filename)
                if isfile(destination)
                    skipped += 1
                    continue
                end
                url = "https://opendap.co-ops.nos.noaa.gov/thredds/fileServer/NOAA/CBOFS/MODELS/$yyyy/$mm/$dd/cbofs.t$(cycle)z.$yyyymmdd.fields.n$hour.nc"
                try
                    println("Downloading $filename")
                    Downloads.download(url, destination)
                    downloaded += 1
                catch err
                    @warn "Download failed" url destination exception=(err, catch_backtrace())
                    isfile(destination) && rm(destination)
                    failed += 1
                end
            end
        end
    end
    println("CBOFS download complete: $downloaded new, $skipped existing, $failed failed.")
end

download_window(
    parse_date_env("START_DATE", DEFAULT_START_DATE),
    parse_date_env("END_DATE", DEFAULT_END_DATE),
    get(ENV, "DATA_DIR", "datafiles"),
)
