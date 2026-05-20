# Alpha Power Regression Checklist

1. Build one subject-level dataset from saved spectra outputs.
- Source: `*_all_electrode_spectra.mat`
- Metric: band-integrated `8-13 Hz` alpha power in dB at `F3`
- Variables: `SubjectID`, `Group`, `AgeYears`, `AlphaPower_dB`, `Electrode`, `SourceFile`

2. Center age across the full analysis sample.
- Compute `age_mean`
- Compute `AgeCentered = AgeYears - age_mean`
- Code group as `GroupCode` with `0 = Control`, `1 = Patient`

3. Fit the primary linear interaction model.
- `AlphaPower_dB ~ AgeCentered * GroupCode`

4. Fit the exploratory quadratic interaction model.
- `AlphaPower_dB ~ AgeCentered * GroupCode + AgeCentered2 * GroupCode`
- Where `AgeCentered2 = AgeCentered.^2`

5. Compute model-comparison metrics manually.
- Residual sum of squares
- AIC
- AICc with `K = p + 1` including residual variance
- Delta AICc

6. Retain Model 1 by default.
- Only promote Model 2 after manual review of diagnostics and biological plausibility
- Use the script output to assess whether Model 2 is a serious candidate

7. Generate residual diagnostics for both models.
- Residuals vs fitted
- Q-Q plot

8. Report coefficient tables and group-specific fitted equations.
- Save coefficient tables for both models
- Save readable control and patient equations for both models

9. Plot fitted curves on the original age axis.
- Raw scatter by group
- Fitted control and patient curves
- 95% confidence bands for fitted means

10. Save all outputs.
- Regression input table
- Model comparison table
- Coefficient tables
- Diagnostics figures
- Fitted-curve figure
- MAT summary file
