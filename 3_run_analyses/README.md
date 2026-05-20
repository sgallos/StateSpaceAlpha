# 3 Run Analyses

This folder contains scripts you run from MATLAB.

Main command:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
RUN_ME_StateSpaceAlpha
```

Key scripts:

- `run_ssm_no_resampling_matrix.m`: main six-run no-subject-resampling
  sensitivity matrix.
- `run_ssm_posterior_no_resampling.m`: primary fixed-hyperparameter posterior
  workflow.
- `run_ss_age_diff_em_fixed_sigma.m`: EM q with fixed empirical sigmaBio.
- `run_ss_age_diff_em_subepoch76.m`: direct 76-row sub-epoch diagnostic.
- `run_ss_age_diff_bootstrap_simultaneous.m`: optional bootstrap sensitivity.

These scripts write full generated outputs to `outputs/`. Curated copies for
review live in `4_figures/` and `5_outputs_data/`.
