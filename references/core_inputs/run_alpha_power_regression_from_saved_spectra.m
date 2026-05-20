%% run_alpha_power_regression_from_saved_spectra
% Regression analysis of alpha power vs age using saved subject-level spectra.
% Uses point-estimate subject data only, not MFDB bootstrap outputs.

%% Step 1: Setup
clear; clc;

scriptRoot = '/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS';
pedsRoot = '/Users/gallo/Documents/peds_cp';
addpath(fullfile(pedsRoot, 'src'));
addpath(genpath(fullfile(pedsRoot, 'src')));

anesthesia = 'sevo';
electrodeLabel = 'F3';
bandName = 'Alpha';
bandRange = [8 13];
ageMin = 10;
ageMax = 15.5;
manualRetainModel2 = false;

outDir = fullfile(scriptRoot, 'outputs', 'alpha_power_regression');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('Electrode: %s\n', electrodeLabel);
fprintf('Band: %s [%g %g] Hz\n', bandName, bandRange(1), bandRange(2));
fprintf('Age range: [%.1f %.1f]\n', ageMin, ageMax);
fprintf('Manual retain Model 2 override: %d\n', manualRetainModel2);

%% Step 2: Build subject-level alpha-power dataset from saved spectra
meta = load_subject_metadata(anesthesia);

controlTable = build_age_table(meta.controls.sortedTable, 'Control');
patientTable = build_age_table(meta.patients.sortedTable, 'Patient');
subjectTable = [controlTable; patientTable];

ageMask = isfinite(subjectTable.AgeYears) & ...
    subjectTable.AgeYears >= ageMin & subjectTable.AgeYears <= ageMax;
subjectTable = sortrows(subjectTable(ageMask, :), {'Group', 'AgeYears'});

subjectTable.SpectraFile = strings(height(subjectTable), 1);
subjectTable.HasSpectraFile = false(height(subjectTable), 1);

for i = 1:height(subjectTable)
    subjectID = subjectTable.SubjectID(i);
    matPath = fullfile(scriptRoot, 'outputs', sprintf('%s_electrode_review', subjectID), ...
        sprintf('%s_all_electrode_spectra.mat', subjectID));
    if isfile(matPath)
        subjectTable.SpectraFile(i) = string(matPath);
        subjectTable.HasSpectraFile(i) = true;
    end
end

missingFiles = subjectTable(~subjectTable.HasSpectraFile, :);
if ~isempty(missingFiles)
    fprintf('\nMissing saved spectra files for:\n');
    disp(missingFiles(:, {'SubjectID', 'Group', 'AgeYears'}));
end

subjectTable = subjectTable(subjectTable.HasSpectraFile, :);

regressionData = table();
regressionData.SubjectID = strings(0, 1);
regressionData.Group = strings(0, 1);
regressionData.AgeYears = zeros(0, 1);
regressionData.AlphaPowerLinear = zeros(0, 1);
regressionData.AlphaPower_dB = zeros(0, 1);
regressionData.Electrode = strings(0, 1);
regressionData.SourceFile = strings(0, 1);

missingElectrode = strings(0, 1);
invalidBand = strings(0, 1);

for i = 1:height(subjectTable)
    subjectID = char(subjectTable.SubjectID(i));
    matPath = char(subjectTable.SpectraFile(i));
    S = load(matPath, 'labels', 'f', 'spectraLinear');

    labels = string(S.labels(:));
    electrodeIdx = find(strcmpi(labels, electrodeLabel), 1, 'first');
    if isempty(electrodeIdx)
        missingElectrode(end + 1, 1) = string(subjectID); %#ok<AGROW>
        continue;
    end

    freq = S.f(:);
    bandMask = freq >= bandRange(1) & freq <= bandRange(2);
    if ~any(bandMask)
        invalidBand(end + 1, 1) = string(subjectID); %#ok<AGROW>
        continue;
    end

    powerLinear = trapz(freq(bandMask), S.spectraLinear(bandMask, electrodeIdx));
    powerLinear = max(powerLinear, eps);
    powerDB = 10 * log10(powerLinear);

    row = table( ...
        string(subjectID), ...
        subjectTable.Group(i), ...
        subjectTable.AgeYears(i), ...
        powerLinear, ...
        powerDB, ...
        string(electrodeLabel), ...
        string(matPath), ...
        'VariableNames', {'SubjectID', 'Group', 'AgeYears', 'AlphaPowerLinear', 'AlphaPower_dB', 'Electrode', 'SourceFile'});

    regressionData = [regressionData; row]; %#ok<AGROW>
end

regressionData = sortrows(regressionData, {'Group', 'AgeYears', 'SubjectID'});

if ~isempty(missingElectrode)
    fprintf('\nMissing electrode %s in:\n', electrodeLabel);
    disp(missingElectrode);
end
if ~isempty(invalidBand)
    fprintf('\nNo valid %s band bins in:\n', bandName);
    disp(invalidBand);
end

fprintf('\nRegression dataset:\n');
disp(regressionData(:, {'SubjectID', 'Group', 'AgeYears', 'AlphaPower_dB'}));

%% Step 3: Center age and code group
nSubjects = height(regressionData);
if nSubjects < 6
    error('Too few subjects for regression: n = %d', nSubjects);
end

ageMean = mean(regressionData.AgeYears, 'omitnan');
regressionData.AgeCentered = regressionData.AgeYears - ageMean;
regressionData.AgeCentered2 = regressionData.AgeCentered .^ 2;
regressionData.GroupCode = double(strcmpi(regressionData.Group, 'Patient'));

fprintf('\nRegression sample size: %d\n', nSubjects);
fprintf('Controls: %d\n', sum(regressionData.GroupCode == 0));
fprintf('Patients: %d\n', sum(regressionData.GroupCode == 1));
fprintf('Age mean used for centering: %.6f years\n', ageMean);

%% Step 4: Fit Model 1 and Model 2
mdl1 = fitlm(regressionData, 'AlphaPower_dB ~ AgeCentered*GroupCode');
mdl2 = fitlm(regressionData, 'AlphaPower_dB ~ AgeCentered*GroupCode + AgeCentered2*GroupCode');

coefTable1 = coefficient_table_with_terms(mdl1);
coefTable2 = coefficient_table_with_terms(mdl2);

%% Step 5: Compute AICc manually and compare models
comparisonTable = [ ...
    build_model_comparison_row(mdl1, 'Model1_Linear'), ...
    build_model_comparison_row(mdl2, 'Model2_Quadratic') ...
    ];
comparisonTable = struct2table(comparisonTable);
comparisonTable.DeltaAICc = comparisonTable.AICc - min(comparisonTable.AICc);
comparisonTable.Retained = false(height(comparisonTable), 1);

model2AICcAdvantage = comparisonTable.AICc(strcmp(comparisonTable.ModelName, 'Model2_Quadratic')) + 2 < ...
    comparisonTable.AICc(strcmp(comparisonTable.ModelName, 'Model1_Linear'));

selectedModelName = "Model1_Linear";
selectionReason = "Model 1 retained by protocol default.";
if manualRetainModel2
    if model2AICcAdvantage
        selectedModelName = "Model2_Quadratic";
        selectionReason = "Manual override applied after diagnostic review; Model 2 also met DeltaAICc > 2 criterion.";
    else
        selectedModelName = "Model1_Linear";
        selectionReason = "Manual override requested, but Model 2 did not meet the DeltaAICc > 2 criterion.";
    end
elseif model2AICcAdvantage
    selectionReason = "Model 2 has DeltaAICc > 2 advantage, but Model 1 remains retained until manual diagnostic review confirms stable, plausible curvature.";
end
comparisonTable.Retained(strcmp(comparisonTable.ModelName, selectedModelName)) = true;

fprintf('\nModel comparison:\n');
disp(comparisonTable);
fprintf('Selection: %s\n', selectedModelName);
fprintf('Reason: %s\n', selectionReason);

%% Step 6: Build readable group-specific equations
model1Equations = build_group_equations_model1(mdl1);
model2Equations = build_group_equations_model2(mdl2);

fprintf('\nModel 1 equations:\n');
fprintf('Controls: %s\n', model1Equations.Controls);
fprintf('Patients: %s\n', model1Equations.Patients);

fprintf('\nModel 2 equations:\n');
fprintf('Controls: %s\n', model2Equations.Controls);
fprintf('Patients: %s\n', model2Equations.Patients);

%% Step 7: Residual diagnostics for both models
figDiag1 = make_diagnostic_figure(mdl1, 'Model 1: Linear Interaction Diagnostics');
figDiag2 = make_diagnostic_figure(mdl2, 'Model 2: Quadratic Interaction Diagnostics');

%% Step 8A: Plot fitted curves with 95%% confidence bands on original age axis
if selectedModelName == "Model2_Quadratic"
    selectedModel = mdl2;
    selectedFormula = 'Quadratic';
else
    selectedModel = mdl1;
    selectedFormula = 'Linear';
end

ageGrid = linspace(min(regressionData.AgeYears), max(regressionData.AgeYears), 200).';
ageGridCentered = ageGrid - ageMean;
controlGrid = make_prediction_table(ageGrid, ageMean, 0);
patientGrid = make_prediction_table(ageGrid, ageMean, 1);

[controlPred, controlCI] = predict(selectedModel, controlGrid);
[patientPred, patientCI] = predict(selectedModel, patientGrid);

ctrlColor = [0.9 0.45 0.1];
patColor = [0.2 0.45 0.85];
ctrlCIColor = [1.0 0.77 0.58];
patCIColor = [0.64 0.80 1.0];

diffColor = [0.1 0.1 0.1];
diffCIColor = [0.55 0.75 0.95];
sigBandColor = [1.0 0.78 0.72];
zeroLineColor = [0.35 0.35 0.35];

figFit = figure('Color', 'w', 'Position', [100 100 1200 560]);
ax = axes(figFit);
hold(ax, 'on');

hCtrlScatter = scatter(ax, regressionData.AgeYears(regressionData.GroupCode == 0), ...
    regressionData.AlphaPower_dB(regressionData.GroupCode == 0), 70, ...
    'MarkerFaceColor', ctrlColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
hPatScatter = scatter(ax, regressionData.AgeYears(regressionData.GroupCode == 1), ...
    regressionData.AlphaPower_dB(regressionData.GroupCode == 1), 70, ...
    'MarkerFaceColor', patColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);

hCtrlBand = fill(ax, [ageGrid; flipud(ageGrid)], [controlCI(:, 1); flipud(controlCI(:, 2))], ctrlCIColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');
hPatBand = fill(ax, [ageGrid; flipud(ageGrid)], [patientCI(:, 1); flipud(patientCI(:, 2))], patCIColor, ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none');

hCtrlLine = plot(ax, ageGrid, controlPred, '-', 'Color', ctrlColor, 'LineWidth', 2.5);
hPatLine = plot(ax, ageGrid, patientPred, '-', 'Color', patColor, 'LineWidth', 2.5);

ax.FontName = 'Helvetica';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.Box = 'off';
grid(ax, 'off');

xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('%s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, sprintf('Alpha Power Regression vs Age (%s retained)', selectedFormula), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hCtrlScatter, hPatScatter, hCtrlBand, hPatBand, hCtrlLine, hPatLine], ...
    {'Controls', 'Patients', 'Controls 95% CI', 'Patients 95% CI', ...
    'Controls fit', 'Patients fit'}, 'Location', 'best');

%% Step 8B: Plot control-only fitted curve
controlEquationText = build_group_equation_annotation(model1Equations, model2Equations, selectedModelName, 'Controls');
controlStatsText = build_group_stats_annotation(selectedModel, 'Controls');
figControl = figure('Color', 'w', 'Position', [100 100 1000 560]);
ax = axes(figControl);
hold(ax, 'on');

hControlScatterOnly = scatter(ax, regressionData.AgeYears(regressionData.GroupCode == 0), ...
    regressionData.AlphaPower_dB(regressionData.GroupCode == 0), 70, ...
    'MarkerFaceColor', ctrlColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
hControlBandOnly = fill(ax, [ageGrid; flipud(ageGrid)], [controlCI(:, 1); flipud(controlCI(:, 2))], ctrlCIColor, ...
    'FaceAlpha', 0.16, 'EdgeColor', 'none');
hControlLineOnly = plot(ax, ageGrid, controlPred, '-', 'Color', ctrlColor, 'LineWidth', 3.0);

ax.FontName = 'Helvetica';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.Box = 'off';
grid(ax, 'off');

xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('%s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, sprintf('Control Alpha Power Regression (%s retained)', selectedFormula), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hControlScatterOnly, hControlLineOnly, hControlBandOnly], ...
    {'Controls', sprintf('%s fit', selectedFormula), '95% CI'}, 'Location', 'best');
text(ax, 0.03, 0.97, controlStatsText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', 'FontSize', 9, ...
    'BackgroundColor', [1 1 1]);
text(ax, 0.03, 0.83, controlEquationText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', 'FontSize', 9, 'Color', [0.35 0.35 0.35], ...
    'BackgroundColor', 'none', 'Interpreter', 'tex');

%% Step 8C: Plot CP-only fitted curve
patientEquationText = build_group_equation_annotation(model1Equations, model2Equations, selectedModelName, 'Patients');
patientStatsText = build_group_stats_annotation(selectedModel, 'Patients');
figPatient = figure('Color', 'w', 'Position', [100 100 1000 560]);
ax = axes(figPatient);
hold(ax, 'on');

hPatientScatterOnly = scatter(ax, regressionData.AgeYears(regressionData.GroupCode == 1), ...
    regressionData.AlphaPower_dB(regressionData.GroupCode == 1), 70, ...
    'MarkerFaceColor', patColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
hPatientBandOnly = fill(ax, [ageGrid; flipud(ageGrid)], [patientCI(:, 1); flipud(patientCI(:, 2))], patCIColor, ...
    'FaceAlpha', 0.16, 'EdgeColor', 'none');
hPatientLineOnly = plot(ax, ageGrid, patientPred, '-', 'Color', patColor, 'LineWidth', 3.0);

ax.FontName = 'Helvetica';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.Box = 'off';
grid(ax, 'off');

xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('%s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, sprintf('CP Alpha Power Regression (%s retained)', selectedFormula), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hPatientScatterOnly, hPatientLineOnly, hPatientBandOnly], ...
    {'Patients', sprintf('%s fit', selectedFormula), '95% CI'}, 'Location', 'best');
text(ax, 0.03, 0.97, patientStatsText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', 'FontSize', 9, ...
    'BackgroundColor', [1 1 1]);
text(ax, 0.03, 0.83, patientEquationText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', 'FontSize', 9, 'Color', [0.35 0.35 0.35], ...
    'BackgroundColor', 'none', 'Interpreter', 'tex');

%% Step 8D: Plot CP-Control difference vs age with 95%% confidence interval
[diffFit, diffCILow, diffCIHigh, diffTerms, diffBeta, diffCov, tCrit, diffPTable] = ...
    compute_difference_curve(selectedModel, ageGridCentered);

sigMask = diffCILow > 0 | diffCIHigh < 0;
sigIntervals = mask_to_intervals(ageGrid, sigMask);
diffEquationText = build_difference_equation_annotation(diffTerms, diffBeta, diffPTable);

figDiff = figure('Color', 'w', 'Position', [100 100 1200 560]);
ax = axes(figDiff);
hold(ax, 'on');

hDiffBand = fill(ax, [ageGrid; flipud(ageGrid)], [diffCILow; flipud(diffCIHigh)], diffCIColor, ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none');
hDiffLine = plot(ax, ageGrid, diffFit, '-', 'Color', diffColor, 'LineWidth', 3.0);
hZero = yline(ax, 0, '--', 'Color', zeroLineColor, 'LineWidth', 1.5);

if any(sigMask)
    yLimits = [min([diffCILow; diffCIHigh]), max([diffCILow; diffCIHigh])];
    if yLimits(1) == yLimits(2)
        yLimits = yLimits + [-1 1];
    end
    hSigBand = plot(ax, nan, nan, 's', 'MarkerFaceColor', sigBandColor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', 'none');
    for i = 1:size(sigIntervals, 1)
        patch(ax, [sigIntervals(i, 1) sigIntervals(i, 2) sigIntervals(i, 2) sigIntervals(i, 1)], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], sigBandColor, ...
            'FaceAlpha', 0.22, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    uistack(hDiffBand, 'top');
    uistack(hDiffLine, 'top');
    uistack(hZero, 'top');
else
    hSigBand = plot(ax, nan, nan, 's', 'MarkerFaceColor', sigBandColor, ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', 'none');
end

ax.FontName = 'Helvetica';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
ax.Box = 'off';
grid(ax, 'off');

xlabel(ax, 'Age (years)');
ylabel(ax, 'Alpha Power Difference (dB)');
title(ax, sprintf('Alpha Power Difference vs Age (%s retained)', selectedFormula), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hSigBand, hDiffBand, hDiffLine, hZero], ...
    {'CI excludes zero', '95% CI', 'CP - Control difference', 'Zero line'}, ...
    'Location', 'best');

text(ax, 0.03, 0.97, diffEquationText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontName', 'Helvetica', 'FontSize', 10, ...
    'BackgroundColor', [1 1 1], 'Interpreter', 'tex');

%% Step 9: Save outputs
inputCsv = fullfile(outDir, 'alpha_power_regression_input.csv');
inputMat = fullfile(outDir, 'alpha_power_regression_input.mat');
comparisonCsv = fullfile(outDir, 'alpha_power_model_comparison.csv');
coef1Csv = fullfile(outDir, 'alpha_power_model1_coefficients.csv');
coef2Csv = fullfile(outDir, 'alpha_power_model2_coefficients.csv');
summaryMat = fullfile(outDir, 'alpha_power_regression_results.mat');

writetable(regressionData, inputCsv);
save(inputMat, 'regressionData', 'ageMean', 'bandName', 'bandRange', 'electrodeLabel', 'ageMin', 'ageMax', 'anesthesia');
writetable(comparisonTable, comparisonCsv);
writetable(coefTable1, coef1Csv);
writetable(coefTable2, coef2Csv);
saveas(figDiag1, fullfile(outDir, 'model1_diagnostics.png'));
saveas(figDiag2, fullfile(outDir, 'model2_diagnostics.png'));
saveas(figFit, fullfile(outDir, 'alpha_power_regression_fit.png'));
saveas(figControl, fullfile(outDir, 'alpha_power_regression_control_fit.png'));
saveas(figPatient, fullfile(outDir, 'alpha_power_regression_cp_fit.png'));
saveas(figDiff, fullfile(outDir, 'alpha_power_difference_vs_age.png'));

save(summaryMat, ...
    'regressionData', 'ageMean', ...
    'mdl1', 'mdl2', ...
    'coefTable1', 'coefTable2', ...
    'comparisonTable', ...
    'selectedModelName', 'selectionReason', ...
    'model1Equations', 'model2Equations', ...
    'controlEquationText', 'patientEquationText', ...
    'ageGrid', 'ageGridCentered', ...
    'controlPred', 'controlCI', 'patientPred', 'patientCI', ...
    'diffFit', 'diffCILow', 'diffCIHigh', 'sigMask', 'sigIntervals', ...
    'diffTerms', 'diffBeta', 'diffCov', 'tCrit', 'diffPTable', 'diffEquationText', ...
    'bandName', 'bandRange', 'electrodeLabel', 'ageMin', 'ageMax', 'anesthesia');

fprintf('\nSaved outputs to %s\n', outDir);

%% Local functions
function T = build_age_table(sortedTable, groupLabel)
    T = table();
    T.SubjectID = string(sortedTable.Properties.RowNames);
    T.Group = repmat(string(groupLabel), height(sortedTable), 1);
    T.AgeYears = sortedTable.('AGE (YEARS)');
end

function out = build_model_comparison_row(mdl, modelName)
    n = mdl.NumObservations;
    p = mdl.NumEstimatedCoefficients;
    K = p + 1;
    rss = sum(mdl.Residuals.Raw .^ 2, 'omitnan');
    sigma2 = rss / n;
    aic = n * log(sigma2) + n * (1 + log(2 * pi)) + 2 * K;
    aicc = aic + (2 * K * (K + 1)) / (n - K - 1);

    out = struct();
    out.ModelName = string(modelName);
    out.NumObservations = n;
    out.NumBetaParameters = p;
    out.K_includingSigma2 = K;
    out.ResidualDF = mdl.DFE;
    out.RSS = rss;
    out.AIC = aic;
    out.AICc = aicc;
    out.AdjustedR2 = mdl.Rsquared.Adjusted;
end

function T = coefficient_table_with_terms(mdl)
    T = mdl.Coefficients;
    T.Term = string(T.Properties.RowNames);
    T = movevars(T, 'Term', 'Before', 1);
end

function eq = build_group_equations_model1(mdl)
    b0 = get_coef(mdl, '(Intercept)');
    bAge = get_coef(mdl, 'AgeCentered');
    bGroup = get_coef(mdl, 'GroupCode');
    bInt = get_coef(mdl, 'AgeCentered:GroupCode');

    eq = struct();
    eq.Controls = build_linear_equation(b0, bAge);
    eq.Patients = build_linear_equation(b0 + bGroup, bAge + bInt);
end

function eq = build_group_equations_model2(mdl)
    b0 = get_coef(mdl, '(Intercept)');
    bAge = get_coef(mdl, 'AgeCentered');
    bGroup = get_coef(mdl, 'GroupCode');
    bInt = get_coef(mdl, 'AgeCentered:GroupCode');
    bAge2 = get_coef(mdl, 'AgeCentered2');
    bAge2Int = get_coef(mdl, 'AgeCentered2:GroupCode');

    eq = struct();
    eq.Controls = build_quadratic_equation(b0, bAge, bAge2);
    eq.Patients = build_quadratic_equation(b0 + bGroup, bAge + bInt, bAge2 + bAge2Int);
end

function value = get_coef(mdl, termName)
    idx = strcmp(mdl.CoefficientNames, termName);
    if any(idx)
        value = mdl.Coefficients.Estimate(idx);
    else
        value = 0;
    end
end

function eqText = build_linear_equation(intercept, slope)
    eqText = sprintf('alpha = %.4f %s %.4f*age_c', intercept, sign_symbol(slope), abs(slope));
end

function eqText = build_quadratic_equation(intercept, slope, quad)
    eqText = sprintf('alpha = %.4f %s %.4f*age_c %s %.4f*age_c^2', ...
        intercept, sign_symbol(slope), abs(slope), sign_symbol(quad), abs(quad));
end

function symbol = sign_symbol(value)
    if value < 0
        symbol = '-';
    else
        symbol = '+';
    end
end

function fig = make_diagnostic_figure(mdl, figTitle)
    fig = figure('Color', 'w', 'Position', [100 100 1150 480]);

    ax1 = subplot(1, 2, 1);
    scatter(ax1, mdl.Fitted, mdl.Residuals.Raw, 60, 'filled', ...
        'MarkerFaceColor', [0.2 0.45 0.85], 'MarkerEdgeColor', 'k');
    hold(ax1, 'on');
    yline(ax1, 0, '--k', 'LineWidth', 1);
    xlabel(ax1, 'Fitted values');
    ylabel(ax1, 'Residuals');
    title(ax1, 'Residuals vs Fitted');
    style_diag_axes(ax1);

    ax2 = subplot(1, 2, 2);
    axes(ax2); %#ok<LAXES>
    qqplot(mdl.Residuals.Raw);
    title('Q-Q Plot');
    style_diag_axes(gca);

    sgtitle(fig, figTitle, 'FontSize', 14, 'FontWeight', 'bold');
end

function style_diag_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.Box = 'off';
    grid(ax, 'off');
end

function predTbl = make_prediction_table(ageGrid, ageMean, groupCode)
    predTbl = table();
    predTbl.AgeYears = ageGrid;
    predTbl.AgeCentered = ageGrid - ageMean;
    predTbl.AgeCentered2 = predTbl.AgeCentered .^ 2;
    predTbl.GroupCode = repmat(groupCode, numel(ageGrid), 1);
end

function [diffFit, diffCILow, diffCIHigh, diffTerms, diffBeta, diffCov, tCrit, diffPTable] = compute_difference_curve(mdl, ageGridCentered)
    if any(strcmp(mdl.CoefficientNames, 'AgeCentered2:GroupCode'))
        diffTerms = {'GroupCode', 'AgeCentered:GroupCode', 'AgeCentered2:GroupCode'};
        X = [ones(numel(ageGridCentered), 1), ageGridCentered(:), ageGridCentered(:) .^ 2];
    else
        diffTerms = {'GroupCode', 'AgeCentered:GroupCode'};
        X = [ones(numel(ageGridCentered), 1), ageGridCentered(:)];
    end

    idx = zeros(1, numel(diffTerms));
    diffPTable = table();
    diffPTable.Term = strings(0, 1);
    diffPTable.PValue = zeros(0, 1);
    for k = 1:numel(diffTerms)
        idx(k) = find(strcmp(mdl.CoefficientNames, diffTerms{k}), 1, 'first');
        if isempty(idx(k))
            error('Required coefficient %s not found in selected model.', diffTerms{k});
        end
        diffPTable = [diffPTable; table(string(diffTerms{k}), mdl.Coefficients.pValue(idx(k)), ...
            'VariableNames', {'Term', 'PValue'})]; %#ok<AGROW>
    end

    diffBeta = mdl.Coefficients.Estimate(idx);
    diffCov = mdl.CoefficientCovariance(idx, idx);
    diffFit = X * diffBeta;
    diffSE = sqrt(diag(X * diffCov * X.'));
    tCrit = tinv(0.975, mdl.DFE);
    diffCILow = diffFit - tCrit * diffSE;
    diffCIHigh = diffFit + tCrit * diffSE;
end

function intervals = mask_to_intervals(xVals, mask)
    intervals = zeros(0, 2);
    if isempty(xVals) || ~any(mask)
        return;
    end

    mask = mask(:);
    xVals = xVals(:);
    dMask = diff([false; mask; false]);
    startIdx = find(dMask == 1);
    endIdx = find(dMask == -1) - 1;

    for i = 1:numel(startIdx)
        intervals(end + 1, :) = [xVals(startIdx(i)), xVals(endIdx(i))]; %#ok<AGROW>
    end
end

function txt = build_group_stats_annotation(mdl, groupName)
    if strcmp(groupName, 'Controls')
        slopeTerm = 'AgeCentered';
    else
        slopeTerm = 'AgeCentered';
        if any(strcmp(mdl.CoefficientNames, 'AgeCentered:GroupCode'))
            slopeTerm = 'AgeCentered:GroupCode';
        end
    end

    coefNames = mdl.CoefficientNames;
    beta = mdl.Coefficients.Estimate;
    covMat = mdl.CoefficientCovariance;
    tCrit = tinv(0.975, mdl.DFE);

    ageIdx = find(strcmp(coefNames, 'AgeCentered'), 1, 'first');
    intIdx = find(strcmp(coefNames, 'AgeCentered:GroupCode'), 1, 'first');

    if strcmp(groupName, 'Controls')
        slope = beta(ageIdx);
        slopeVar = covMat(ageIdx, ageIdx);
    else
        slope = beta(ageIdx);
        slopeVar = covMat(ageIdx, ageIdx);
        if ~isempty(intIdx)
            slope = slope + beta(intIdx);
            slopeVar = slopeVar + covMat(intIdx, intIdx) + 2 * covMat(ageIdx, intIdx);
        end
    end

    slopeSe = sqrt(max(slopeVar, 0));
    slopeLow = slope - tCrit * slopeSe;
    slopeHigh = slope + tCrit * slopeSe;

    fitLabel = 'Linear fit';
    if any(strcmp(coefNames, 'AgeCentered2'))
        fitLabel = 'Quadratic fit';
    end

    txt = sprintf('%s\nBeta_age = %.2f dB/year\n95%% CI [%.2f, %.2f]\nR^2 = %.2f', ...
        fitLabel, slope, slopeLow, slopeHigh, mdl.Rsquared.Ordinary);
end

function txt = build_group_equation_annotation(model1Equations, model2Equations, selectedModelName, groupName)
    if selectedModelName == "Model2_Quadratic"
        if strcmp(groupName, 'Controls')
            txt = model2Equations.Controls;
        else
            txt = model2Equations.Patients;
        end
    else
        if strcmp(groupName, 'Controls')
            txt = model1Equations.Controls;
        else
            txt = model1Equations.Patients;
        end
    end
end

function txt = build_difference_equation_annotation(diffTerms, diffBeta, diffPTable)
    if numel(diffTerms) == 2
        offset = diffBeta(1);
        slope = diffBeta(2);
        pGroup = diffPTable.PValue(strcmp(diffPTable.Term, 'GroupCode'));
        txt = sprintf('\\Delta\\alpha = %.2f %s %.2f*age_c\\newline p(group) = %.3f', ...
            offset, sign_symbol(slope), abs(slope), pGroup);
    else
        offset = diffBeta(1);
        slope = diffBeta(2);
        quad = diffBeta(3);
        pGroup = diffPTable.PValue(strcmp(diffPTable.Term, 'GroupCode'));
        txt = sprintf('\\Delta\\alpha = %.2f %s %.2f*age_c %s %.2f*age_c^2\\newline p(group) = %.3f', ...
            offset, sign_symbol(slope), abs(slope), sign_symbol(quad), abs(quad), pGroup);
    end
end
