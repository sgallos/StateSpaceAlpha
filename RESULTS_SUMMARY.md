# Results Summary

## Headline Finding

The no-subject-resampling SSM analysis supports a lower frontal alpha-power
trajectory in CP subjects than in controls during sevoflurane anesthesia.

Primary fixed-hyperparameter run:

- Mean delta(a): -4.85 dB
- Delta(a) range: -5.98 to -2.74 dB
- Age range where the 95% posterior credible band excludes zero: 12.20 to
  14.87 years

Here `delta(a)` means CP minus Control alpha power. Negative values mean CP
has lower alpha power than Control.

See:

- [4_figures/subject_level_alpha_vs_age.png](4_figures/subject_level_alpha_vs_age.png)
- [4_figures/primary_delta_trajectory.png](4_figures/primary_delta_trajectory.png)
- [4_figures/sensitivity_all_six_overlay.png](4_figures/sensitivity_all_six_overlay.png)
- [5_outputs_data/ssm_no_resampling_matrix_summary.csv](5_outputs_data/ssm_no_resampling_matrix_summary.csv)

## Six-Run No-Resampling Matrix

All runs use the SSM posterior covariance for inference. Subjects are not
resampled. The runs differ only in how the hyperparameters are chosen.

| Strategy name | What it does | q_f0 | q_delta | sigmaBio | Converged | Boundary | Zero excluded age range |
| --- | --- | ---: | ---: | ---: | --- | --- | --- |
| Fixed-Heuristic | Fixed qInit, fixed empirical sigmaBio | 2.97 | 2.97 | 14.51 | yes | no | 12.20-14.87 yr |
| EM-Smoothness | EM q, fixed empirical sigmaBio | 14.04 | 25.54 | 14.51 | no | no | 13.93-14.87 yr |
| EM-Joint | Joint EM, collapsed 19-row table | 15.26 | 36.52 | 15.75 | no | no | 14.03-14.58 yr |
| EM-Joint-Subepoch | Joint EM, direct 76-row sub-epoch table | 1e6 | 1e6 | 1.32 | yes | yes | diagnostic failure |
| EM-Biological | Fixed qInit, EM sigmaBio | 2.97 | 2.97 | 15.72 | yes | no | 12.20-14.87 yr |
| Profile-Likelihood | Profile-likelihood q, fixed empirical sigmaBio | 3.86 | 0.01 | 14.51 | yes | yes | 12.20-15.20 yr |

Interpretation:

- Fixed-Heuristic, EM-Biological, and Profile-Likelihood are the cleanest
  no-resampling sensitivity runs.
- EM-Smoothness and EM-Joint are useful diagnostics but did not converge
  cleanly.
- EM-Joint-Subepoch should not be used as a primary result because both
  smoothness parameters hit the ceiling.
- Five of the six runs support the same qualitative conclusion: CP alpha
  power is lower than Control alpha power, with separation strongest in the
  older part of the age range.

## Main Figures

- [subepoch_consistency_scatter.png](4_figures/subepoch_consistency_scatter.png):
  four 30 s alpha estimates per subject, used to check recording stability.
- [alpha_table_subepoch_vs_120s_comparison.png](4_figures/alpha_table_subepoch_vs_120s_comparison.png):
  compares the original 120 s table to the collapsed 4 x 30 s table.
- [subject_level_alpha_vs_age.png](4_figures/subject_level_alpha_vs_age.png):
  raw 19-subject alpha-vs-age scatter before modeling.
- [primary_delta_trajectory.png](4_figures/primary_delta_trajectory.png):
  primary fixed-hyperparameter delta(a) posterior mean and credible band.
- [subjects_with_trajectories.png](4_figures/subjects_with_trajectories.png):
  Control and CP posterior trajectories with the 19 subject-level alpha values
  and their within-subject/sub-epoch error bars overlaid.
- [sensitivity_all_six_overlay.png](4_figures/sensitivity_all_six_overlay.png):
  posterior means from all six no-resampling runs.
- [sensitivity_all_six_bands.png](4_figures/sensitivity_all_six_bands.png):
  each run shown separately with its posterior credible band.
- [sensitivity_profile_loglik_heatmap.png](4_figures/sensitivity_profile_loglik_heatmap.png):
  Profile-Likelihood surface for `q_f0` and `q_delta`.
- [em_convergence_diagnostics.png](4_figures/em_convergence_diagnostics.png):
  EM iteration histories explaining which EM strategies converged cleanly and
  which should be treated as diagnostics.

## Scientific Reading

The SSM result is consistent with a CP-related reduction in frontal alpha
power under sevoflurane. The exact age range of statistical separation is
sensitive to hyperparameter choice, which is expected at n = 19. The robust
claim is not "the effect starts at exactly age 12.20"; the robust claim is
"the CP-Control difference is negative across the age range and is most
clearly separated between roughly 12 and 15 years."
