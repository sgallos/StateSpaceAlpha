# 5 Outputs Data

Curated small CSV outputs for review.

- `subepoch_consistency_summary.csv`: per-subject QC summary for the four
  30 s alpha estimates.
- `alpha_table_subepoch_vs_120s_comparison.csv`: compares the collapsed
  sub-epoch alpha table against the original 120 s MFDB table.
- `ssm_no_resampling_matrix_summary.csv`: one row per no-resampling strategy,
  including hyperparameters, convergence flags, boundary flags, and the age
  range where the credible band excludes zero.
- `ssm_no_resampling_matrix_trajectories.csv`: trajectory-level output for
  all six no-resampling runs.
- `ssm_posterior_no_resampling_primary_trajectory.csv`: trajectory output for
  the primary fixed-hyperparameter run.
- `ssm_posterior_no_resampling_sensitivity_summary.csv`: fixed-grid
  sensitivity summary from the primary no-resampling workflow.
- `ssm_posterior_no_resampling_sigmaBio_subject_diagnostics.csv`: subject-level
  diagnostics used to estimate empirical `sigmaBio`.

Large `.mat` files and full generated outputs remain in `outputs/`, which is
ignored by git.
