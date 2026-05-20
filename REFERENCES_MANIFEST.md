# References Manifest

This repository is intentionally a skeleton workspace. It contains the reference files we currently have locally, plus scaffolds for the new state-space analysis. It does not need to contain every historical source file from every prior project.

## Included Core Inputs

These were copied from the current `MATLAB_Multitaper_Hz_Domain_BTS` pipeline:

- `references/core_inputs/mfdb_config.m`
- `references/core_inputs/run_subject_mfdb_bootstrap_for_id.m`
- `references/core_inputs/mfdb_subject_registry.csv`
- `references/core_inputs/CS29_F3_mfdb_fixture.mat`
- `references/core_inputs/plot_mfdb_alpha_power_vs_age.m`
- `references/core_inputs/plot_mfdb_alpha_power_difference_vs_age.m`
- `references/core_inputs/run_alpha_power_regression_from_saved_spectra.m`

## Included Documentation

- `references/docs/MFDB_CHECKLIST.md`
- `references/docs/ALPHA_POWER_REGRESSION_CHECKLIST.md`

## Included SSM Templates

These are reference-only templates, not production code for this project:

- `references/ssm_templates/ay_fft_filter.m`
- `references/ssm_templates/EM_parameters_kim_2018.m` copied from `/Users/gallo/projects/SSMT/EM_parameters.m`
- `references/ssm_templates/computeAMI_SMI.m`
- `references/ssm_templates/kim_ssmt_2018/EM_parameters.m`
- `references/ssm_templates/kim_ssmt_2018/SS_MT.m`
- `references/ssm_templates/kim_ssmt_2018/SS_ST.m`
- `references/ssm_templates/kim_ssmt_2018/main.m`
- `references/ssm_templates/kim_ssmt_2018/multitaper.m`
- `references/ssm_templates/kim_ssmt_2018/periodogram.m`
- `references/ssm_templates/kim_ssmt_2018/README.md`
- `references/ssm_templates/kim_ssmt_2018/SED10.mat`

## New Scaffolds

These are the files to implement in this repo:

- `SS_age_diff.m`
- `run_ss_age_diff.m`
- `smoke_test_alpha_extraction.m`
- `hgam_sensitivity.R`
- `SIMULATION_PLAN.md`

## Not Included Unless Added Later

The repo does not currently include:

- all 19 full per-subject MFDB `.mat` files
- the full 134 MB `CS29_F3_mfdb.mat`; only the small fixture is included
- the full `peds_cp` source tree
- the full `MATLAB_Multitaper_Hz_Domain_BTS` analysis tree
- any external paper PDFs
- full Kim SSMT git history; a small source snapshot is included under `references/ssm_templates/kim_ssmt_2018`

That is intentional. The purpose is to keep the workspace small and readable while preserving enough reference material to implement and review the state-space alpha-age analysis.
