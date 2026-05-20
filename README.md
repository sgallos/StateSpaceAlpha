# StateSpaceAlpha

State-space model for frontal alpha power in children with cerebral palsy
(CP) and neurotypical controls during sevoflurane anesthesia.

## What This Project Does

This project analyzes frontal alpha power at F3 in children aged 10 to
15.5 years. The cohort has 11 CP subjects and 8 controls. The scientific
question is:

Does CP change the alpha-power age trajectory under sevoflurane?

The analysis uses a varying-coefficient state-space model. The key output is
delta(a), the CP minus Control alpha-power difference as a function of age.
Credible bands come directly from the smoother posterior covariance. No
subject resampling is used in the primary no-resampling analysis.

## Main Result

Across six no-subject-resampling SSM runs, five runs support the same
qualitative result: CP subjects have lower frontal alpha power than controls,
by about 5 dB on average. The posterior credible band excludes zero somewhere
in the 12 to 15 year age range in all five non-degenerate runs.

The direct 76-row sub-epoch EM diagnostic also excludes zero, but its
smoothness parameters hit the ceiling. Treat that run as a diagnostic failure,
not as a primary result.

Read:

- [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) for the result table and figure guide.
- [METHODS.md](METHODS.md) for the model in plain English.
- [SIMULATION_PLAN.md](SIMULATION_PLAN.md) for validation work still needed.

## How To Run The Main Analysis

Open MATLAB:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
RUN_ME_StateSpaceAlpha
```

The wrapper has clearly labeled sections and switches at the top. By default,
it runs the primary no-resampling posterior workflow, the six-run
no-resampling matrix, and refreshes curated review files.

The key curated copies are stored in:

- `4_figures/`
- `5_outputs_data/`

## Folder Guide

- `1_data_preparation/`: turns MFDB outputs into alpha-power tables.
- `2_state_space_model/`: contains the Kalman filter, RTS smoother, EM code,
  and variance-estimation helper.
- `3_run_analyses/`: scripts that run the analyses and produce outputs.
- `4_figures/`: curated figures for review.
- `5_outputs_data/`: curated small CSV outputs for review.
- `references/`: core input references and analysis notes from the original
  multitaper pipeline.
- `reference_code/`: external or template SSM code kept for comparison.
- `reference_papers/`: paper pointers and methodological anchors.

## Recommended Reading Order

1. `RESULTS_SUMMARY.md`
2. `METHODS.md`
3. `4_figures/README.md`
4. `RUN_ME_StateSpaceAlpha.m`
5. `3_run_analyses/run_ssm_no_resampling_matrix.m`
6. `2_state_space_model/SS_age_diff.m`

## Notes

The full generated `outputs/` folder is intentionally ignored by git. The
curated figures and small CSV summaries are tracked separately so a
collaborator can inspect the result without rerunning the full analysis.
