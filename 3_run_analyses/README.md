# 3 Analysis Scripts

This folder contains scripts you run from MATLAB.

Main command:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
RUN_ME_StateSpaceAlpha
```

Key scripts:

- `run_power_ssm_scalar_r.m`: power-band SSM workflow with one EM-fit total
  observation-noise variance `r`. The current default runs alpha (8-12 Hz)
  on the collapsed sub-epoch table, compares four q treatments, and writes
  trajectory, group-difference, q-version, EM-diagnostic, and Profile-q
  heatmap figures.
- `run_approach1_mfdb_noise_only.m`: MFDB observation-noise-only analysis.
  Runs Path A (120 s MFDB variance) and Path B (sub-epoch SEM variance) with
  both likelihood-selected grid q and fixed heuristic q.
- `plot_approach1_path_comparison.m`: four-panel comparison of Approach 1
  Path A/Path B and grid-q/fixed-q fits.
- `run_ssm_no_resampling_matrix.m`: main six-run no-subject-resampling
  sensitivity matrix.
- `plot_em_convergence_diagnostics.m`: convergence diagnostic for the four
  EM-based strategies in the sensitivity matrix.
- `plot_subject_level_alpha_vs_age.m`: raw 19-subject alpha-vs-age scatter
  before modeling.
- `plot_subjects_with_trajectories.m`: subject-level alpha values with
  Control and CP posterior trajectories overlaid.
- `run_ssm_posterior_no_resampling.m`: primary fixed-hyperparameter posterior
  workflow.
- `run_ss_age_diff_em_fixed_sigma.m`: EM q with fixed empirical sigmaBio.
- `run_ss_age_diff_em_subepoch76.m`: direct 76-row sub-epoch diagnostic.
- `run_ss_age_diff_bootstrap_simultaneous.m`: optional bootstrap sensitivity.

These scripts write full generated outputs to `outputs/`. Curated copies for
review live in `4_figures/` and `5_outputs_data/`.
