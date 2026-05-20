%% run_ssm_posterior_no_resampling
% Primary SSM posterior workflow with no subject resampling.
%
% This script treats the SSM smoother posterior as the inferential output.
% Hyperparameters are fixed, not EM-estimated. The CP-Control difference
% trajectory delta(a) and its 95% credible band come directly from:
%
%   smoothedState(3, k) +/- 1.96 * sqrt(smoothedCovariance(3, 3, k))
%
% The sensitivity analysis changes fixed hyperparameters, not subjects.

clear; clc;

%% Step 1: User-facing settings
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

collapsedAlphaTableFileName = 'alpha_table_for_ssm_subepochs_collapsed.csv';
subepochAlphaTableFileName = 'alpha_table_for_ssm_subepochs.csv';

% Fixed sigmaBio estimation from the 76-row sub-epoch table.
% Options: "trend_residual_sem", "trend_residual_mfdb", "between_only".
sigmaBioEstimator = "trend_residual_sem";

% Options: "constant", "linear", "quadratic".
sigmaBioTrendModel = "linear";

% Primary process-noise scale. qInit is computed from the alpha range and
% age span, then multiplied by these sensitivity values.
primaryQScaleF0 = 1;
primaryQScaleDelta = 1;
qScaleGrid = [0.1 1 10];

% Sensitivity grid for fixed biological variance. The empirical estimate is
% inserted below after sigmaBio is computed.
extraSigmaBioSensitivityValues = [0.5 1 2 5 10];

%% Step 2: Load alpha tables
collapsedAlphaTableFile = fullfile(outDir, collapsedAlphaTableFileName);
if ~isfile(collapsedAlphaTableFile)
    error(['Missing collapsed alpha table:\n%s\n\n' ...
        'Run collapse_subepoch_alpha_table_for_ssm first.'], collapsedAlphaTableFile);
end

subepochAlphaTableFile = fullfile(outDir, subepochAlphaTableFileName);
upstreamSubepochAlphaTableFile = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', subepochAlphaTableFileName);
if ~isfile(subepochAlphaTableFile)
    if isfile(upstreamSubepochAlphaTableFile)
        copyfile(upstreamSubepochAlphaTableFile, subepochAlphaTableFile);
        fprintf('Copied sub-epoch alpha table into StateSpaceAlpha outputs:\n  %s\n', ...
            subepochAlphaTableFile);
    else
        error(['Missing sub-epoch alpha table:\n%s\n\n' ...
            'Run smoke_test_alpha_extraction_subepochs first.'], subepochAlphaTableFile);
    end
end

alphaTable = readtable(collapsedAlphaTableFile, 'TextType', 'string');
subepochTable = readtable(subepochAlphaTableFile, 'TextType', 'string');

fprintf('Loaded collapsed alpha table: %s\n', collapsedAlphaTableFile);
fprintf('Loaded sub-epoch alpha table: %s\n', subepochAlphaTableFile);
fprintf('Collapsed rows: %d | Sub-epoch rows: %d\n', height(alphaTable), height(subepochTable));

%% Step 3: Fixed hyperparameters
[sigmaBioEstimate, sigmaBioDiagnostics] = estimate_sigmaBio_from_subepochs(subepochTable, ...
    'Estimator', sigmaBioEstimator, ...
    'TrendModel', sigmaBioTrendModel, ...
    'SigmaBioFloor', 0.01, ...
    'Verbose', true);

qInit = compute_q_init(alphaTable.ageYears, alphaTable.alpha_dB);
primaryQF0 = primaryQScaleF0 * qInit;
primaryQDelta = primaryQScaleDelta * qInit;

fprintf('\nPrimary fixed hyperparameters\n');
fprintf('  qInit: %.6g\n', qInit);
fprintf('  q_f0: %.6g\n', primaryQF0);
fprintf('  q_delta: %.6g\n', primaryQDelta);
fprintf('  sigmaBio: %.6g dB^2\n', sigmaBioEstimate);

%% Step 4: Primary smoother fit
primaryResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', primaryQF0, ...
    'processNoiseQDelta', primaryQDelta, ...
    'biologicalVariance', sigmaBioEstimate, ...
    'verbose', true);

primaryTrajectory = primaryResult.trajectory;
deltaMean = primaryTrajectory.deltaMean_dB;
deltaSD = primaryTrajectory.deltaSD_dB;
deltaCILow = primaryTrajectory.deltaCILow_dB;
deltaCIHigh = primaryTrajectory.deltaCIHigh_dB;
ageYears = primaryTrajectory.ageYears;

bandExcludesZero = any(deltaCILow > 0) || any(deltaCIHigh < 0);
positiveAges = ageYears(deltaCILow > 0);
negativeAges = ageYears(deltaCIHigh < 0);

fprintf('\nPrimary delta(a) posterior summary\n');
fprintf('  Mean delta across ages: %.3f dB\n', mean(deltaMean, 'omitnan'));
fprintf('  Delta range: %.3f to %.3f dB\n', min(deltaMean), max(deltaMean));
fprintf('  Delta 95%% CI low range: %.3f to %.3f dB\n', min(deltaCILow), max(deltaCILow));
fprintf('  Delta 95%% CI high range: %.3f to %.3f dB\n', min(deltaCIHigh), max(deltaCIHigh));
fprintf('  Band excludes zero at any age: %d\n', bandExcludesZero);
if ~isempty(positiveAges)
    fprintf('  Positive exclusion ages: %.3f to %.3f years\n', min(positiveAges), max(positiveAges));
end
if ~isempty(negativeAges)
    fprintf('  Negative exclusion ages: %.3f to %.3f years\n', min(negativeAges), max(negativeAges));
end

%% Step 5: Sensitivity grid over fixed q and sigmaBio
sigmaBioGrid = unique([extraSigmaBioSensitivityValues sigmaBioEstimate 2 * sigmaBioEstimate]);
nQ = numel(qScaleGrid);
nSigma = numel(sigmaBioGrid);
sensitivityResults = cell(nQ, nSigma);
sensitivitySummary = table('Size', [nQ * nSigma, 9], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double','logical'}, ...
    'VariableNames', {'qScale','sigmaBio','qF0','qDelta','meanDelta', ...
    'minDelta','maxDelta','maxBandWidth','bandExcludesZero'});

rowIdx = 0;
for i = 1:nQ
    for j = 1:nSigma
        rowIdx = rowIdx + 1;
        thisQF0 = qScaleGrid(i) * qInit;
        thisQDelta = qScaleGrid(i) * qInit;
        thisSigmaBio = sigmaBioGrid(j);

        sensitivityResults{i, j} = SS_age_diff(alphaTable, ...
            'processNoiseQF0', thisQF0, ...
            'processNoiseQDelta', thisQDelta, ...
            'biologicalVariance', thisSigmaBio, ...
            'verbose', false);

        thisTrajectory = sensitivityResults{i, j}.trajectory;
        thisBandWidth = thisTrajectory.deltaCIHigh_dB - thisTrajectory.deltaCILow_dB;

        sensitivitySummary.qScale(rowIdx) = qScaleGrid(i);
        sensitivitySummary.sigmaBio(rowIdx) = thisSigmaBio;
        sensitivitySummary.qF0(rowIdx) = thisQF0;
        sensitivitySummary.qDelta(rowIdx) = thisQDelta;
        sensitivitySummary.meanDelta(rowIdx) = mean(thisTrajectory.deltaMean_dB, 'omitnan');
        sensitivitySummary.minDelta(rowIdx) = min(thisTrajectory.deltaMean_dB);
        sensitivitySummary.maxDelta(rowIdx) = max(thisTrajectory.deltaMean_dB);
        sensitivitySummary.maxBandWidth(rowIdx) = max(thisBandWidth);
        sensitivitySummary.bandExcludesZero(rowIdx) = ...
            any(thisTrajectory.deltaCILow_dB > 0) || any(thisTrajectory.deltaCIHigh_dB < 0);
    end
end

%% Step 6: Save outputs
primaryTrajectoryFile = fullfile(outDir, 'ssm_posterior_no_resampling_primary_trajectory.csv');
sensitivitySummaryFile = fullfile(outDir, 'ssm_posterior_no_resampling_sensitivity_summary.csv');
sigmaBioSubjectDiagnosticsFile = fullfile(outDir, 'ssm_posterior_no_resampling_sigmaBio_subject_diagnostics.csv');
matFile = fullfile(outDir, 'ssm_posterior_no_resampling_results.mat');

writetable(primaryTrajectory, primaryTrajectoryFile);
writetable(sensitivitySummary, sensitivitySummaryFile);
writetable(sigmaBioDiagnostics.subjectTable, sigmaBioSubjectDiagnosticsFile);
save(matFile, 'primaryResult', 'sensitivityResults', 'sensitivitySummary', ...
    'sigmaBioEstimate', 'sigmaBioDiagnostics', 'sigmaBioEstimator', ...
    'sigmaBioTrendModel', 'qInit', ...
    'primaryQF0', 'primaryQDelta', 'qScaleGrid', 'sigmaBioGrid', ...
    'primaryTrajectory', 'alphaTable', 'subepochTable');

fprintf('\nSaved primary trajectory to %s\n', primaryTrajectoryFile);
fprintf('Saved sensitivity summary to %s\n', sensitivitySummaryFile);
fprintf('Saved sigmaBio subject diagnostics to %s\n', sigmaBioSubjectDiagnosticsFile);
fprintf('Saved MAT results to %s\n', matFile);

%% Step 7: Primary posterior delta plot
figPrimary = figure('Color', 'w', 'Position', [100 100 950 520]);
ax = axes(figPrimary);
hold(ax, 'on');

hBand = plot_ci_band(ax, ageYears, deltaCILow, deltaCIHigh, [0.70 0.80 0.90], 0.60);
hDelta = plot(ax, ageYears, deltaMean, '-k', 'LineWidth', 2.5);
hZero = yline(ax, 0, '--', 'Color', [0.40 0.40 0.40], 'LineWidth', 1.2);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control delta(a) (dB)');
title(ax, sprintf('SSM posterior delta(a), no resampling: q = %.3g, sigmaBio = %.3g', ...
    qInit, sigmaBioEstimate));
legend(ax, [hBand, hDelta, hZero], ...
    {'95% posterior credible band', 'Posterior mean', 'Zero'}, 'Location', 'best');
style_ssm_axes(ax);

primaryPlotFile = fullfile(outDir, 'ssm_posterior_delta_primary_no_resampling.png');
saveas(figPrimary, primaryPlotFile);
fprintf('Saved primary posterior plot to %s\n', primaryPlotFile);

%% Step 8: Sensitivity overlay plot
figSensitivity = figure('Color', 'w', 'Position', [100 100 1120 560]);
ax = axes(figSensitivity);
hold(ax, 'on');
yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);

colors = lines(nQ * nSigma);
plotIdx = 0;
for i = 1:nQ
    for j = 1:nSigma
        plotIdx = plotIdx + 1;
        thisTrajectory = sensitivityResults{i, j}.trajectory;
        plot(ax, thisTrajectory.ageYears, thisTrajectory.deltaMean_dB, ...
            '-', 'LineWidth', 1.4, 'Color', colors(plotIdx, :), ...
            'DisplayName', sprintf('q x %.1f, sigmaBio %.3g', qScaleGrid(i), sigmaBioGrid(j)));
    end
end

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control delta(a) (dB)');
title(ax, 'SSM posterior delta(a) sensitivity, no resampling');
legend(ax, 'Location', 'eastoutside');
style_ssm_axes(ax);

sensitivityPlotFile = fullfile(outDir, 'ssm_posterior_delta_sensitivity_no_resampling.png');
saveas(figSensitivity, sensitivityPlotFile);
fprintf('Saved sensitivity plot to %s\n', sensitivityPlotFile);

fprintf('\nNo-resampling posterior workflow complete.\n');

%% Local helpers
function qInit = compute_q_init(ageYears, alphaDB)
    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange ^ 2) / max(ageSpan ^ 3, eps);
end

function hBand = plot_ci_band(ax, xVals, ciLow, ciHigh, bandColor, faceAlpha)
    xVals = xVals(:);
    ciLow = ciLow(:);
    ciHigh = ciHigh(:);
    validMask = isfinite(xVals) & isfinite(ciLow) & isfinite(ciHigh);

    hBand = patch(ax, ...
        [xVals(validMask); flipud(xVals(validMask))], ...
        [ciLow(validMask); flipud(ciHigh(validMask))], ...
        bandColor, ...
        'EdgeColor', 'none', 'FaceAlpha', faceAlpha);
    uistack(hBand, 'bottom');
end

function style_ssm_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
