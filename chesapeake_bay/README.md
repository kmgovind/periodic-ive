# Chesapeake Bay simulation study

This directory contains the ACC 2027 Chesapeake Bay case study. The
release-reproduction workflow is `cb_sigma_comparison.jl`; it runs the
adaptive IVE planner and a matched constant-speed baseline for
\(\sigma_\zeta \in \{2,3,4\}\) PSU.

## Data

The simulations use NOAA Chesapeake Bay Operational Forecast System (CBOFS)
surface-salinity fields. The paper window is 16--21 February 2026 UTC. Inputs
must be placed in `datafiles/` and named:

```text
chesapeake_salinity_YYYYMMDD_tCCz_nFFF.nc
```

Here `CC` is the model cycle and `FFF` is the three-digit forecast hour. The
data directory is ignored by Git. Download the archive with:

```bash
cd chesapeake_bay
julia --project=.. download_cbofs_data.jl
```

The script is restart-safe and queries NOAA's daily THREDDS catalog. The full
window is roughly 7 GB in the currently archived data. Retain NOAA attribution
in any redistributed archive.

## Reproduce the paper experiment

From this directory:

```bash
julia --project=.. cb_sigma_comparison.jl
```

This produces, under `sigma_comparison/`:

- one JLD2 result archive per adaptive and constant-speed run;
- `sigma_comparison_summary.csv` with clarity and speed metrics;
- `speed_vs_full_path_by_sigma.pdf`; and
- `speed_vs_time_by_strategy.pdf`.

The default configuration matches the submitted
manuscript: 200 path nodes, 60 s simulation step, \(S_0=0.005\), sensing
length scale \(\sigma=1500\) m, clarity-decay rate \(\alpha=0.001\),
nominal speed 1.75 m/s, and mean available power 750 W. There is no
maximum-lap-time constraint. `VehicleParams` retains a 0.01 m/s numerical
guard so the travel-time parameterization remains finite.

The script is standalone: the canonical figures remain in
`sigma_comparison/`. To additionally copy them to a separate manuscript
checkout, set `PAPER_FIGURE_DIR` to that checkout's figure directory.

To regenerate figures and CSV from an existing result archive without
rerunning the costly simulations:

```bash
PLOT_ONLY=true julia --project=.. cb_sigma_comparison.jl
```

## Other scripts

- `cb_sim.jl` is the production dynamic simulator retained from development.
  Do not modify it merely to reproduce the submitted figures; use
  `cb_sigma_comparison.jl` instead.
- `cb_max_lap_time_sweep.jl` is a separate, optional sweep of a maximum lap
  duration in 24-hour increments.
- `importance_weight_decay_sensitivity.jl` is an earlier screening study that
  explicitly imposes a 0.5 m/s speed floor. It is useful for historical
  sensitivity analysis but is not the final paper configuration.
- `cb_data.ipynb`, `path_plan.ipynb`, and the other notebooks are exploratory
  provenance. The scripts above are the maintained reproduction interfaces.

## Result interpretation

The total-clarity values are comparable between IVE and constant speed within
one \(\sigma_\zeta\) row. They should not be compared directly across
\(\sigma_\zeta\) values, because changing this parameter changes the target
weighting objective itself.
