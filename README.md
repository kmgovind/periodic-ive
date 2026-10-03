# Information Value of Energy: reproducibility code

This repository contains the Julia implementation and numerical studies for
the ACC 2027 submission, *Optimally Trading Energy for Information through
Speed Trajectory Optimization of Sustainably Powered Mobile Robots*.

The release-ready workflows are intentionally script-based. The older
notebooks are retained as exploratory records, but are not required to
reproduce the paper figures or numerical results.

## Setup

Use Julia with the pinned environment in this directory:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The Chesapeake Bay experiments require the NOAA CBOFS NetCDF files in
`chesapeake_bay/datafiles/`. They are not version-controlled because the
five-day archive is about 7 GB. See
[`chesapeake_bay/README.md`](chesapeake_bay/README.md) for an exact download
command, the expected filename convention, and the full reproduction path.

## Paper workflows

| Goal | Script | Output |
| --- | --- | --- |
| Illustrative-example grid study | `simple_sensitivity/toy_problem.jl` | `simple_sensitivity/grid_independence_sweep.pdf` |
| Chesapeake Bay paper experiment | `chesapeake_bay/cb_sigma_comparison.jl` | results and figures under `chesapeake_bay/sigma_comparison/` |
| Regenerate the two Chesapeake Bay speed figures from archived results | `PLOT_ONLY=true julia --project=.. cb_sigma_comparison.jl` (run in `chesapeake_bay/`) | `sigma_comparison/*.pdf` |
| Optional maximum-lap-time sensitivity | `chesapeake_bay/cb_max_lap_time_sweep.jl` | `chesapeake_bay/max_lap_time_sweep/` |
| Optional target-weight sensitivity with a 0.5 m/s floor | `chesapeake_bay/importance_weight_decay_sensitivity.jl` | timestamped local output files |

The checked-in `Manifest.toml` records the environment used for the study.

## Repository layout

- `src/`: clarity, energy, CBOFS-data, and optimization utilities.
- `simple_sensitivity/`: the illustrative numerical validation.
- `chesapeake_bay/`: data download, dynamic simulations, figure regeneration,
  and sensitivity studies.
- `verification/`: derivations and exploratory verification material; not part
  of the paper reproduction workflow.
- `gulf_stream/` and `habs/`: exploratory work outside the ACC paper scope.

## Data and result release policy

The repository ignores NetCDF inputs, JLD2 result archives, generated figures,
and videos. For a public release, publish the exact CBOFS input archive and
the `sigma_comparison/` result archive in a DOI-backed repository (for example
Zenodo), then add its DOI to `chesapeake_bay/README.md`. This keeps the Git
repository compact while still allowing exact figure regeneration.

## License and citation

Add a license file and a `CITATION.cff` before making the repository public.
The current repository does not yet declare a redistribution license.
