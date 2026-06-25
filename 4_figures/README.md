# Figures

Read in this order. Each figure includes its key takeaway.

## subepoch_consistency_scatter.png

For each subject, the 4 sub-epoch alpha values are plotted side by side.
Takeaway: most subjects show tight clustering, meaning the 120 s epochs were
recorded under reasonably stable conditions. A few subjects may show larger
within-subject spread and should be inspected if the sensitivity result
depends on them.

## alpha_table_subepoch_vs_120s_comparison.png

Compares `mfdb_var` per subject under two approaches: the original 120 s MFDB
bootstrap vs the empirical variance of 4 sub-epochs.
Takeaway: the sub-epoch approach gives larger and more realistic observation
variances because it reflects within-recording physiology, not just multitaper
estimation noise.

## subject_level_alpha_vs_age.png

The raw data the model sees: 19 subjects, colored by group.
Takeaway: there is no obvious smooth age trend within either group. CP values
are shifted lower than Control values, and within-group scatter is substantial
relative to any age structure.

## primary_delta_trajectory.png

The headline result. CP-minus-Control alpha trajectory across age. The shaded
region is the 95% posterior credible band.
Takeaway: the trajectory is consistently negative (CP < Control) across the
whole age range, with the band excluding zero between about ages 12.2 and
14.9 years in the primary fixed-hyperparameter run.

## subjects_with_trajectories.png

The 19 subjects plotted as points, overlaid against the posterior Control
trajectory `f0(a)` and CP trajectory `f0(a) + delta(a)`.
Takeaway: this is the concrete data-plus-model view. The Control trajectory is
visibly higher than the CP trajectory, while the subject points show which
children sit above or below their group's posterior expectation.

## sensitivity_all_six_overlay.png

Six hyperparameter selection strategies overlaid on the same axes.
Takeaway: the trajectories agree qualitatively. The CP < Control finding is
robust to how smoothness and biological variance are chosen.

## sensitivity_all_six_bands.png

Same six strategies, now with credible bands instead of just curves.
Takeaway: band width varies, but the qualitative finding is preserved. The
age range where the band excludes zero shifts slightly across strategies.

## sensitivity_profile_loglik_heatmap.png

Log-likelihood as a function of `(q_f0, q_delta)` from the
Profile-Likelihood strategy.
Takeaway: the surface favors very smooth `delta(a)` and indicates that
smoothness is not sharply determined by the data at n = 19, which is expected
for this cohort size.

## em_convergence_diagnostics.png

Iteration histories for the four EM-based strategies. Each column is one
strategy, with log-likelihood, smoothness parameters, and `sigmaBio` shown
from top to bottom.
Takeaway: this figure shows why some EM strategies are useful diagnostics
rather than primary analyses: clean runs should have monotonic log-likelihood,
settled `q` values, and no boundary hit.

## approach1_four_panel_comparison.png

Approach 1 sets observation noise exactly equal to `mfdb_var` and sets
`sigma2_bio = 0`. The four panels compare Path A (120 s MFDB variance) and
Path B (sub-epoch SEM variance), each under likelihood-selected grid `q` and
fixed heuristic `q`.
Takeaway: the grid-q fits are more flexible and can land on the q-grid
boundary, while the fixed-q fits impose smoother trajectories. This is the
direct check of what the SSM says when MFDB variance is the only observation
noise.

## approach1_PathA_120s_*.png and approach1_PathB_subepoch_*.png

Path-specific Approach 1 outputs. Each path has a posterior delta-band plot,
a fixed-q delta-band plot, and a log-likelihood heatmap over the q grid.
Path A uses pure 120 s MFDB spectral-estimation variance. Path B uses the
sub-epoch-derived SEM variance from the collapsed 4 x 30 s table.
