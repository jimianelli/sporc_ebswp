# EBS pollock SPoRC assessment

This repository is a starter assessment wrapper for EBS walleye pollock using
the local `SPoRC` R package. It mirrors the nearby `fims` and `scam` assessment
repos: configuration is kept in `config/`, reusable R code in `R/`, data
conditioning in `data-raw/`, model runs in `analysis/`, and reporting in
`reporting/`.

## Layout

- `config/base.yml`: base SPoRC model and path configuration.
- `config/scenarios.yml`: scenario grid for selectivity random effects.
- `R/build_inputs.R`: converts conditioned pollock data into SPoRC inputs.
- `R/extract_results.R`: extracts SPoRC time series for comparison.
- `data-raw/build_ebs_pollock_data.R`: builds `data/ebs_pollock.rds` from the
  2024 ADMB data file.
- `analysis/run_scenarios.R`: fits configured SPoRC scenarios.
- `analysis/compare_admb.R`: compares the SPoRC base fit with
  `/Users/jim/_mymods/pollock/results/derived/admb_2024.rds`.

## Workflow

```sh
Rscript data-raw/build_ebs_pollock_data.R
Rscript analysis/run_scenarios.R
Rscript analysis/compare_admb.R
```

The comparison script writes:

- `analysis/comparisons/sporc_vs_admb_2024.rds`
- `analysis/comparisons/sporc_vs_admb_2024.csv`

The default ADMB source root is `/Users/jim/_mymods/pollock/admb`. Override it
for data rebuilds with `EBS_POLLOCK_ADMB=/path/to/admb`.
