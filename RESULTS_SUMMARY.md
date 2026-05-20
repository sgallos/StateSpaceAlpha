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

- [4_figures/primary_delta_trajectory.png](4_figures/primary_delta_trajectory.png)
- [4_figures/no_resampling_matrix_overlay.png](4_figures/no_resampling_matrix_overlay.png)
- [5_outputs_data/ssm_no_resampling_matrix_summary.csv](5_outputs_data/ssm_no_resampling_matrix_summary.csv)

## Six-Run No-Resampling Matrix

All runs use the SSM posterior covariance for inference. Subjects are not
resampled. The runs differ only in how the hyperparameters are chosen.

| Run | Strategy | q_f0 | q_delta | sigmaBio | Converged | Boundary | Zero excluded age range |
| --- | --- | ---: | ---: | ---: | --- | --- | --- |
| A | Fixed qInit, fixed empirical sigmaBio | 2.97 | 2.97 | 14.51 | yes | no | 12.20-14.87 yr |
| B | EM q, fixed empirical sigmaBio | 14.04 | 25.54 | 14.51 | no | no | 13.93-14.87 yr |
| C | Joint EM, collapsed 19-row table | 15.26 | 36.52 | 15.75 | no | no | 14.03-14.58 yr |
| D | Joint EM, direct 76-row sub-epoch table | 1e6 | 1e6 | 1.32 | yes | yes | diagnostic failure |
| E | Fixed qInit, EM sigmaBio | 2.97 | 2.97 | 15.72 | yes | no | 12.20-14.87 yr |
| F | Profile-likelihood q, fixed empirical sigmaBio | 3.86 | 0.01 | 14.51 | yes | yes | 12.20-15.20 yr |

Interpretation:

- A, E, and F are the cleanest no-resampling sensitivity runs.
- B and C are useful diagnostics but did not converge cleanly.
- D should not be used as a primary result because both smoothness parameters
  hit the ceiling.
- Five of the six runs support the same qualitative conclusion: CP alpha
  power is lower than Control alpha power, with separation strongest in the
  older part of the age range.

## Main Figures

- [primary_delta_trajectory.png](4_figures/primary_delta_trajectory.png):
  primary fixed-hyperparameter delta(a) posterior mean and credible band.
- [no_resampling_matrix_overlay.png](4_figures/no_resampling_matrix_overlay.png):
  posterior means from all six no-resampling runs.
- [no_resampling_matrix_delta_bands.png](4_figures/no_resampling_matrix_delta_bands.png):
  each run shown separately with its posterior credible band.
- [profile_likelihood_heatmap.png](4_figures/profile_likelihood_heatmap.png):
  Run F profile-likelihood surface for `q_f0` and `q_delta`.

## Scientific Reading

The SSM result is consistent with a CP-related reduction in frontal alpha
power under sevoflurane. The exact age range of statistical separation is
sensitive to hyperparameter choice, which is expected at n = 19. The robust
claim is not "the effect starts at exactly age 12.20"; the robust claim is
"the CP-Control difference is negative across the age range and is most
clearly separated between roughly 12 and 15 years."
