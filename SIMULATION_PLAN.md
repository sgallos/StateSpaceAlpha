# Simulation Plan

The state-space model should not be trusted on the clinical cohort until it passes basic calibration and identifiability checks. This file is the checklist to run before making scientific claims from `SS_age_diff.m`.

## 1. Input Recovery Check

Goal: verify that the new analysis reads the existing MFDB outputs correctly.

- Load `mfdb_subject_registry.csv`.
- Load the example `CS29_F3_mfdb.mat`.
- Recompute alpha power using the project convention: integrate 8-13 Hz in linear power, then convert to dB.
- Confirm the value matches the existing alpha-power script output for the same subject.

Pass criterion: field names, dimensions, and alpha-power calculation match the reference pipeline.

## 2. Known Linear Difference

Goal: test whether the SSM can recover a simple CP-control difference that varies linearly with age.

- Simulate control alpha power as a linear age function plus noise.
- Simulate CP alpha power with a known offset and known slope difference.
- Fit the IWP SSM.
- Check whether the posterior/smoothed difference covers the true difference curve.

Pass criterion: the estimated CP-control difference follows the true curve and uncertainty bands have sensible width.

## 3. Known Smooth Nonlinear Difference

Goal: test the benefit of the smoothness-prior model.

- Simulate a smooth nonlinear age difference, for example a broad sigmoid or quadratic-like departure.
- Fit the same SSM without changing hyperparameters by hand.
- Compare against the OLS baseline and HGAM sensitivity model.

Pass criterion: SSM recovers the broad trajectory without inventing oscillatory structure.

## 4. Null Difference / Type I Behavior

Goal: test false positive behavior.

- Simulate CP and control subjects from the same age trajectory.
- Keep the real cohort ages and group sizes.
- Repeat many times.
- Record how often the fitted difference excludes zero over broad age regions.

Pass criterion: the method should not repeatedly find structured differences under the null.

## 5. Small-N Stress Test

Goal: follow the Auger-Methe warning that simple SSMs can fail or become weakly identifiable.

- Repeat simulations with the real group sizes.
- Repeat with one or two subjects removed.
- Repeat with uneven age support.
- Track convergence, parameter estimates, and uncertainty inflation.

Pass criterion: the model reports uncertainty when data support is weak instead of producing overconfident trajectories.

## 6. Parametric Bootstrap

Goal: assess model calibration after fitting the real data.

- Fit the model to the observed cohort.
- Simulate replicated datasets from the fitted model.
- Refit the model to each replicate.
- Compare observed residuals, trajectory smoothness, and age-difference features to the replicated distribution.

Pass criterion: observed data should look plausible under the fitted model. If not, revise the model before interpreting the CP-control difference.

## 7. Sensitivity Checks

Run at least these comparisons:

- Different process-noise starting values.
- Wider and narrower process-noise bounds.
- Mean vs median MFDB-derived subject summaries if both are carried forward.
- HGAM sensitivity model in `hgam_sensitivity.R`.

The final write-up should report whether the CP-control alpha-age trajectory is stable across these choices.

