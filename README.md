# StateSpaceAlpha

This repository is a clean workspace for developing a state-space model of alpha-band power trajectories in pediatric cerebral palsy (CP) and control subjects under sevoflurane anesthesia. The working cohort comes from the existing `MATLAB_Multitaper_Hz_Domain_BTS` pipeline and contains children in the 10 to 15.5 year age range, with per-subject F3 multitaper frequency-domain bootstrap (MFDB) spectra generated from one 120 s artifact-free maintenance epoch.

The scientific question is whether CP modifies the age trajectory of frontal alpha power. The existing OLS regression code provides a useful baseline, but it imposes a parametric age trend on a small cohort and does not directly use the MFDB uncertainty pools. The goal here is to build an age-varying state-space model that estimates the CP-control alpha-power difference as a smooth function of age while preserving the uncertainty structure inherited from the spectral analysis.

The methodological extension is a varying-coefficient integrated Wiener process (IWP) state-space model for the age axis, with HGAM treated as a sensitivity analysis. The state-space model is intended to be the primary smoothness-prior analysis, while the MFDB windowed figures remain a transparent nonparametric companion. This is defensible if the model is treated as a structured smoothness assumption, not as proof created by smoothing; it needs simulation, identifiability checks, and sensitivity comparisons before being used for primary claims.

## Layout

```text
StateSpaceAlpha/
  README.md
  SIMULATION_PLAN.md
  smoke_test_alpha_extraction.m
  run_ssm_residual_sanity_check.m
  SS_age_pooled_iwp.m
  SS_age_diff.m
  SS_age_diff_em.m
  SS_age_diff_bootstrap_simultaneous.m
  run_ss_age_diff.m
  run_ss_age_diff_em.m
  run_ss_age_diff_bootstrap_simultaneous.m
  hgam_sensitivity.R
  references/
    core_inputs/
      mfdb_config.m
      run_subject_mfdb_bootstrap_for_id.m
      mfdb_subject_registry.csv
      CS29_F3_mfdb_fixture.mat
      plot_mfdb_alpha_power_vs_age.m
      plot_mfdb_alpha_power_difference_vs_age.m
      run_alpha_power_regression_from_saved_spectra.m
    ssm_templates/
      ay_fft_filter.m
      EM_parameters_kim_2018.m
      computeAMI_SMI.m
      kim_ssmt_2018/
        EM_parameters.m
        SS_MT.m
        SS_ST.m
        multitaper.m
        periodogram.m
        main.m
        SED10.mat
        README.md
    docs/
      MFDB_CHECKLIST.md
      ALPHA_POWER_REGRESSION_CHECKLIST.md
  outputs/
```

## Reference Files

- `references/core_inputs/mfdb_config.m`: locked MFDB constants, including F3, 120 s epoch, TW, K, bootstrap count, and frequency passband.
- `references/core_inputs/run_subject_mfdb_bootstrap_for_id.m`: exact subject-level MFDB generation and saved field names.
- `references/core_inputs/mfdb_subject_registry.csv`: source table for subject IDs, group labels, ages, and per-subject MFDB file paths.
- `references/core_inputs/CS29_F3_mfdb_fixture.mat`: small CS29 fixture with metadata, alpha summary, dimensions, and a few sample values. The full 134 MB subject `.mat` is intentionally not committed.
- `references/core_inputs/plot_mfdb_alpha_power_vs_age.m`: current alpha-power age-window plotting conventions.
- `references/core_inputs/plot_mfdb_alpha_power_difference_vs_age.m`: current alpha-power difference plotting conventions.
- `references/core_inputs/run_alpha_power_regression_from_saved_spectra.m`: current OLS baseline that the SSM analysis is meant to replace or contextualize.
- `references/ssm_templates/ay_fft_filter.m`: main matrix-form Kalman filter, RTS smoother, and EM template.
- `references/ssm_templates/EM_parameters_kim_2018.m`: Kim SSMT EM reference implementation copied from the cloned `sgallos/SSMT` repository.
- `references/ssm_templates/computeAMI_SMI.m`: simpler one-dimensional smoother reference.
- `references/ssm_templates/kim_ssmt_2018/`: small copied snapshot of the cloned `sgallos/SSMT` repository for local reference.

## Paper Pointers

- Kim, Behr, Ba, Brown. State-space multitaper time-frequency analysis. PNAS, 2018. DOI: `10.1073/pnas.1702877115`.
- Kim, Ba, Brown. A multitaper frequency-domain bootstrap method. IEEE Signal Processing Letters, 2018. DOI: `10.1109/LSP.2018.2876606`.
- Song, Chakravarty, Brown. A smoother state space multitaper spectrogram. EMBC, 2018. DOI: `10.1109/EMBC.2018.8512190`.
- Auger-Methe et al. State-space models' dirty little secrets: even simple linear Gaussian models can have estimation problems. Scientific Reports, 2016. DOI: `10.1038/srep26677`.
- Auger-Methe et al. A guide to state-space modeling of ecological time series. Ecological Monographs, 2021. DOI: `10.1002/ecm.1470`.
- Cornelissen et al. Electroencephalographic markers of brain development during sevoflurane anaesthesia in children up to 3 years old. British Journal of Anaesthesia, 2018. DOI: `10.1016/j.bja.2018.01.037`.

## Intended Workflow

1. Run `smoke_test_alpha_extraction.m` to confirm alpha extraction matches the existing MFDB pipeline.
2. Read `mfdb_subject_registry.csv` and load each subject's MFDB `.mat` file.
3. Compute alpha power as already done in the existing MFDB code: integrate 8-13 Hz in linear power, then convert to dB.
4. Build the response vector as subject-level alpha power in dB, with group and age as covariates.
5. Use `SS_age_pooled_iwp.m` as the preserved Step 2 pooled 2-D IWP filter and RTS smoother reference.
6. Run `run_ss_age_diff.m` for the Step 3 4-D baseline + CP-control difference model.
7. Run `run_ssm_residual_sanity_check.m` to compare fitted residual scatter against MFDB observation uncertainty.
8. Run `run_ss_age_diff_em.m` for the Step 4 EM-estimated IWP smoothness model, including the additive biological variance term.
9. Run `run_ss_age_diff_bootstrap_simultaneous.m` for Step 5 two-stage subject + MFDB bootstrap simultaneous confidence bands.
10. Validate the model with parametric bootstrap and identifiability checks before interpreting the fitted trajectory.
11. Run `hgam_sensitivity.R` as a separate smoothness-prior cross-check.

Step 3 run command:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
run_ss_age_diff
```

Current `run_ss_age_diff.m` runs the Step 3 model and writes:

- `outputs/ssm_step3_age_difference_results.mat`
- `outputs/ssm_step3_age_difference_trajectory.csv`
- `outputs/ssm_step3_baseline_cp_overlay.png`
- `outputs/ssm_step3_delta_primary.png`
- `outputs/ssm_step3_delta_q_sensitivity.png`

Step 4 run command:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
run_ss_age_diff_em
```

Current `run_ss_age_diff_em.m` runs the EM-estimated 4-D model and writes:

- `outputs/ssm_step4_em_age_difference_results.mat`
- `outputs/ssm_step4_em_age_difference_trajectory.csv`
- `outputs/ssm_step4_em_multistart_summary.csv`
- `outputs/ssm_step4_em_delta_primary.png`
- `outputs/ssm_step4_em_loglik_history.png`
- `outputs/ssm_step4_em_hyperparameter_history.png`
- `outputs/ssm_step4_em_vs_step3_sensitivity.png`, if Step 3 outputs are present

Step 5 run command:

```matlab
cd('/Users/gallo/projects/StateSpaceAlpha')
run_ss_age_diff_bootstrap_simultaneous
```

Current `run_ss_age_diff_bootstrap_simultaneous.m` runs the two-stage subject-within-group + MFDB-row bootstrap with fixed q and biological-variance values and writes:

- `outputs/boot_alpha_dB_cache.mat`
- `outputs/ssm_step5_bootstrap_simultaneous_results.mat`
- `outputs/ssm_step5_bootstrap_band_summary.csv`
- `outputs/ssm_step5_bootstrap_delta_distribution.csv`
- `outputs/ssm_step5_delta_simultaneous_band.png`
- `outputs/ssm_step5_pointwise_vs_simultaneous_band.png`
- `outputs/ssm_step5_bootstrap_trajectory_diagnostic.png`
