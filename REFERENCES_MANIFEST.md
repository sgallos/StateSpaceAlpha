# References Manifest

This repository is organized for human review. It contains the reference
inputs we currently have locally, the state-space analysis code, curated
figures, and small summary CSVs. It does not include the full historical
source trees or full per-subject MFDB `.mat` outputs.

## Analysis Code

Data preparation:

- `1_data_preparation/smoke_test_alpha_extraction.m`
- `1_data_preparation/collapse_subepoch_alpha_table_for_ssm.m`

State-space model:

- `2_state_space_model/SS_age_diff.m`
- `2_state_space_model/SS_age_diff_em.m`
- `2_state_space_model/SS_age_diff_bootstrap_simultaneous.m`
- `2_state_space_model/SS_age_pooled_iwp.m`
- `2_state_space_model/estimate_sigmaBio_from_subepochs.m`

Analysis runners:

- `RUN_ME_StateSpaceAlpha.m`
- `3_run_analyses/run_ss_age_diff.m`
- `3_run_analyses/run_ssm_posterior_no_resampling.m`
- `3_run_analyses/run_ssm_no_resampling_matrix.m`
- `3_run_analyses/run_ss_age_diff_em.m`
- `3_run_analyses/run_ss_age_diff_em_fixed_sigma.m`
- `3_run_analyses/run_ss_age_diff_em_subepoch76.m`
- `3_run_analyses/run_ss_age_diff_bootstrap_simultaneous.m`
- `3_run_analyses/run_ssm_residual_sanity_check.m`
- `3_run_analyses/hgam_sensitivity.R`

## Included Core Inputs

These were copied from the current `MATLAB_Multitaper_Hz_Domain_BTS`
pipeline:

- `references/core_inputs/mfdb_config.m`
- `references/core_inputs/run_subject_mfdb_bootstrap_for_id.m`
- `references/core_inputs/mfdb_subject_registry.csv`
- `references/core_inputs/CS29_F3_mfdb_fixture.mat`
- `references/core_inputs/plot_mfdb_alpha_power_vs_age.m`
- `references/core_inputs/plot_mfdb_alpha_power_difference_vs_age.m`
- `references/core_inputs/run_alpha_power_regression_from_saved_spectra.m`
- `references/core_inputs/subepoch_mfdb_pipeline/run_subject_mfdb_bootstrap_subepochs_for_id.m`
- `references/core_inputs/subepoch_mfdb_pipeline/run_subject_mfdb_bootstrap_subepochs_batch.m`
- `references/core_inputs/subepoch_mfdb_pipeline/smoke_test_alpha_extraction_subepochs.m`
- `references/core_inputs/subepoch_mfdb_pipeline/check_subepoch_consistency.m`

## Included Documentation

- `METHODS.md`
- `RESULTS_SUMMARY.md`
- `SIMULATION_PLAN.md`
- `references/docs/MFDB_CHECKLIST.md`
- `references/docs/ALPHA_POWER_REGRESSION_CHECKLIST.md`

## Curated Review Outputs

Figures:

- `4_figures/primary_delta_trajectory.png`
- `4_figures/no_resampling_matrix_overlay.png`
- `4_figures/no_resampling_matrix_delta_bands.png`
- `4_figures/profile_likelihood_heatmap.png`
- `4_figures/em_loglik_history.png`

Small CSV outputs:

- `5_outputs_data/ssm_no_resampling_matrix_summary.csv`
- `5_outputs_data/ssm_no_resampling_matrix_trajectories.csv`
- `5_outputs_data/ssm_posterior_no_resampling_primary_trajectory.csv`
- `5_outputs_data/ssm_posterior_no_resampling_sensitivity_summary.csv`
- `5_outputs_data/ssm_posterior_no_resampling_sigmaBio_subject_diagnostics.csv`

## Included Reference Code

These are reference-only templates, not the production analysis:

- `reference_code/ay_fft_filter.m`
- `reference_code/EM_parameters_kim_2018.m`
- `reference_code/computeAMI_SMI.m`
- `reference_code/kim_ssmt_2018/EM_parameters.m`
- `reference_code/kim_ssmt_2018/SS_MT.m`
- `reference_code/kim_ssmt_2018/SS_ST.m`
- `reference_code/kim_ssmt_2018/main.m`
- `reference_code/kim_ssmt_2018/multitaper.m`
- `reference_code/kim_ssmt_2018/periodogram.m`
- `reference_code/kim_ssmt_2018/README.md`
- `reference_code/kim_ssmt_2018/SED10.mat`

## Not Included Unless Added Later

- all 19 full per-subject MFDB `.mat` files
- the full 134 MB `CS29_F3_mfdb.mat`; only the small fixture is included
- the full `peds_cp` source tree
- the full `MATLAB_Multitaper_Hz_Domain_BTS` analysis tree
- external paper PDFs
- full Kim SSMT git history; only a small source snapshot is included

That is intentional. The goal is a compact, readable workspace that preserves
enough context to review and rerun the state-space alpha-age analysis.
