%% run_ss_age_diff_em_fixed_sigma
% Step 4 Path B: fix sigmaBio from sub-epoch structure and estimate only q.
%
% This path avoids joint EM estimation of sigmaBio and the IWP smoothness
% parameters. The 76-row sub-epoch table is used to estimate a fixed
% between-subject variance component, then the 19-row collapsed table is fit
% with sigmaBio held fixed while EM estimates q_f0 and q_delta.

clear; clc;

%% Step 1: Locate input tables
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

subepochTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs.csv');
upstreamSubepochTableFile = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'alpha_table_for_ssm_subepochs.csv');
if ~isfile(subepochTableFile)
    if isfile(upstreamSubepochTableFile)
        copyfile(upstreamSubepochTableFile, subepochTableFile);
        fprintf('Copied sub-epoch table into StateSpaceAlpha outputs:\n  %s\n', subepochTableFile);
    else
        error('Missing sub-epoch table. Run smoke_test_alpha_extraction_subepochs first.');
    end
end

collapsedTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs_collapsed.csv');
if ~isfile(collapsedTableFile)
    error(['Missing collapsed sub-epoch table:\n%s\n\n' ...
        'Run collapse_subepoch_alpha_table_for_ssm first.'], collapsedTableFile);
end

subepochTable = readtable(subepochTableFile, 'TextType', 'string');
alphaTable = readtable(collapsedTableFile, 'TextType', 'string');

fprintf('Loaded sub-epoch table: %s\n', subepochTableFile);
fprintf('Loaded collapsed SSM table: %s\n', collapsedTableFile);

%% Step 2: Estimate fixed sigmaBio from the 76-row table
sigmaBioEstimator = "trend_residual_sem";
trendModel = "linear";
[sigmaBioFixed, sigmaDiagnostics] = estimate_sigmaBio_from_subepochs(subepochTable, ...
    'Estimator', sigmaBioEstimator, ...
    'TrendModel', trendModel, ...
    'SigmaBioFloor', 0.01, ...
    'Verbose', true);

sigmaDiagnosticsTable = sigmaDiagnostics.subjectTable;
sigmaDiagnosticsTableFile = fullfile(outDir, 'ssm_fixed_sigma_subepoch_subject_diagnostics.csv');
writetable(sigmaDiagnosticsTable, sigmaDiagnosticsTableFile);
fprintf('Saved sigmaBio subject diagnostics to %s\n', sigmaDiagnosticsTableFile);

%% Step 3: Primary EM fit with fixed sigmaBio
emOptions = struct();
emOptions.maxIter = 100;
emOptions.tolerance = 1e-5;
emOptions.qInitMultipliers = [0.01 0.1 1 10 100];
emOptions.useApproximateMstep = false;

fprintf('\n=== Primary fixed-sigma run: sigmaBio = %.6g dB^2 ===\n', sigmaBioFixed);
emResultsPrimary = SS_age_diff_em(alphaTable, ...
    'fixSigmaBio', true, ...
    'initialSigmaBio', sigmaBioFixed, ...
    'maxIter', emOptions.maxIter, ...
    'tolerance', emOptions.tolerance, ...
    'qInitMultipliers', emOptions.qInitMultipliers, ...
    'useApproximateMstep', emOptions.useApproximateMstep, ...
    'verbose', true);

primaryTrajectory = emResultsPrimary.trajectory;
primarySummary = build_multistart_summary(emResultsPrimary);

%% Step 4: Sensitivity grid over fixed sigmaBio values
sigmaBioGrid = unique([0.5 1 2 sigmaBioFixed 5 10]);
emResultsSensitivity = cell(numel(sigmaBioGrid), 1);
sensitivitySummary = table('Size', [numel(sigmaBioGrid), 9], ...
    'VariableTypes', {'double','double','double','double','double','double','logical','logical','logical'}, ...
    'VariableNames', {'sigmaBioFixed','qF0','qDelta','finalLogLikelihood','emIterations', ...
    'bestStartIndex','converged','logLikelihoodIsMonotone','hitQBoundary'});

for i = 1:numel(sigmaBioGrid)
    fprintf('\n=== Sensitivity fixed-sigma run %d/%d: sigmaBio = %.6g dB^2 ===\n', ...
        i, numel(sigmaBioGrid), sigmaBioGrid(i));
    emResultsSensitivity{i} = SS_age_diff_em(alphaTable, ...
        'fixSigmaBio', true, ...
        'initialSigmaBio', sigmaBioGrid(i), ...
        'maxIter', emOptions.maxIter, ...
        'tolerance', emOptions.tolerance, ...
        'qInitMultipliers', [0.1 1 10], ...
        'useApproximateMstep', emOptions.useApproximateMstep, ...
        'verbose', false);

    bestStart = emResultsSensitivity{i}.bestStart;
    sensitivitySummary.sigmaBioFixed(i) = sigmaBioGrid(i);
    sensitivitySummary.qF0(i) = emResultsSensitivity{i}.q_f0_em;
    sensitivitySummary.qDelta(i) = emResultsSensitivity{i}.q_delta_em;
    sensitivitySummary.finalLogLikelihood(i) = emResultsSensitivity{i}.finalLogLikelihood;
    sensitivitySummary.emIterations(i) = bestStart.emIterations;
    sensitivitySummary.bestStartIndex(i) = emResultsSensitivity{i}.bestStartIndex;
    sensitivitySummary.converged(i) = emResultsSensitivity{i}.converged;
    sensitivitySummary.logLikelihoodIsMonotone(i) = bestStart.logLikelihoodIsMonotone;
    sensitivitySummary.hitQBoundary(i) = emResultsSensitivity{i}.hitQFloor || emResultsSensitivity{i}.hitQCeiling;
end

%% Step 5: Save tables and MAT results
primaryTrajectoryFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_primary_trajectory.csv');
primarySummaryFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_primary_multistart_summary.csv');
sensitivitySummaryFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_sensitivity_summary.csv');
matFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_results.mat');

writetable(primaryTrajectory, primaryTrajectoryFile);
writetable(primarySummary, primarySummaryFile);
writetable(sensitivitySummary, sensitivitySummaryFile);
save(matFile, 'alphaTable', 'subepochTable', 'sigmaBioFixed', 'sigmaDiagnostics', ...
    'trendModel', 'emResultsPrimary', 'sigmaBioGrid', 'emResultsSensitivity', ...
    'sensitivitySummary', 'emOptions', 'subepochTableFile', 'collapsedTableFile');

fprintf('\nSaved fixed-sigma primary trajectory to %s\n', primaryTrajectoryFile);
fprintf('Saved fixed-sigma primary multistart summary to %s\n', primarySummaryFile);
fprintf('Saved fixed-sigma sensitivity summary to %s\n', sensitivitySummaryFile);
fprintf('Saved fixed-sigma MAT results to %s\n', matFile);

%% Step 6: Primary delta figure
figPrimary = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figPrimary);
hold(ax, 'on');

hBand = plot_ci_band(ax, primaryTrajectory.ageYears, ...
    primaryTrajectory.deltaCILow_dB, primaryTrajectory.deltaCIHigh_dB, ...
    [0.78 0.80 0.86], 0.42);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hDelta = plot(ax, primaryTrajectory.ageYears, primaryTrajectory.deltaMean_dB, ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, sprintf('Fixed sigmaBio EM: CP-Control alpha-power difference (sigmaBio %.3g)', ...
    sigmaBioFixed));
legend(ax, [hBand, hDelta, hZero], {'Delta 95% CI', 'Delta mean', 'Zero line'}, ...
    'Location', 'best');
style_ssm_axes(ax);

primaryPlotFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_delta_primary.png');
saveas(figPrimary, primaryPlotFile);
fprintf('Saved fixed-sigma primary delta plot to %s\n', primaryPlotFile);

%% Step 7: Sensitivity overlay figure
figSensitivity = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figSensitivity);
hold(ax, 'on');
yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);

colors = lines(numel(sigmaBioGrid));
for i = 1:numel(sigmaBioGrid)
    thisTrajectory = emResultsSensitivity{i}.trajectory;
    plot(ax, thisTrajectory.ageYears, thisTrajectory.deltaMean_dB, ...
        '-', 'Color', colors(i, :), 'LineWidth', 1.8, ...
        'DisplayName', sprintf('sigmaBio %.3g', sigmaBioGrid(i)));
end

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, 'Fixed sigmaBio sensitivity: delta(a)');
legend(ax, 'Location', 'eastoutside');
style_ssm_axes(ax);

sensitivityPlotFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_delta_sensitivity.png');
saveas(figSensitivity, sensitivityPlotFile);
fprintf('Saved fixed-sigma sensitivity plot to %s\n', sensitivityPlotFile);

%% Step 8: Primary convergence figures
figLogLik = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figLogLik);
hold(ax, 'on');
for i = 1:numel(emResultsPrimary.startResults)
    thisStart = emResultsPrimary.startResults{i};
    plot(ax, 1:numel(thisStart.logLikHistory), thisStart.logLikHistory, ...
        '-o', 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'DisplayName', sprintf('start %.2g', emResultsPrimary.qInitMultipliers(i)));
end
xlabel(ax, 'EM iteration');
ylabel(ax, 'Log-likelihood');
title(ax, 'Fixed sigmaBio EM: primary log-likelihood histories');
legend(ax, 'Location', 'best');
style_ssm_axes(ax);

logLikPlotFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_loglik_history.png');
saveas(figLogLik, logLikPlotFile);
fprintf('Saved fixed-sigma log-likelihood plot to %s\n', logLikPlotFile);

figQ = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figQ);
hold(ax, 'on');
for i = 1:numel(emResultsPrimary.startResults)
    thisStart = emResultsPrimary.startResults{i};
    iterVals = 1:size(thisStart.qHistory, 1);
    semilogy(ax, iterVals, thisStart.qHistory(:, 1), '-', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q f0 start %.2g', emResultsPrimary.qInitMultipliers(i)));
    semilogy(ax, iterVals, thisStart.qHistory(:, 2), '--', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q delta start %.2g', emResultsPrimary.qInitMultipliers(i)));
end
xlabel(ax, 'EM iteration');
ylabel(ax, 'q value');
title(ax, 'Fixed sigmaBio EM: primary q histories');
legend(ax, 'Location', 'eastoutside');
set(ax, 'YScale', 'log');
style_ssm_axes(ax);

qPlotFile = fullfile(outDir, 'ssm_step4_em_fixed_sigma_q_history.png');
saveas(figQ, qPlotFile);
fprintf('Saved fixed-sigma q-history plot to %s\n', qPlotFile);

%% Step 9: Console summary
fprintf('\nFixed sigmaBio Step 4 summary\n');
fprintf('  trend model: %s\n', trendModel);
fprintf('  sigmaBio fixed: %.6g dB^2\n', sigmaBioFixed);
fprintf('  q_f0_em: %.6g\n', emResultsPrimary.q_f0_em);
fprintf('  q_delta_em: %.6g\n', emResultsPrimary.q_delta_em);
fprintf('  final log-likelihood: %.6f\n', emResultsPrimary.finalLogLikelihood);
fprintf('  converged: %d\n', emResultsPrimary.converged);
fprintf('  best start log-likelihood monotone: %d\n', ...
    emResultsPrimary.bestStart.logLikelihoodIsMonotone);
disp(sensitivitySummary);

%% Local helpers
function summaryTable = build_multistart_summary(emResults)
    nStarts = numel(emResults.startResults);
    startMultiplier = emResults.qInitMultipliers(:);
    initialQF0 = nan(nStarts, 1);
    initialQDelta = nan(nStarts, 1);
    initialSigmaBio = nan(nStarts, 1);
    finalQF0 = nan(nStarts, 1);
    finalQDelta = nan(nStarts, 1);
    finalSigmaBio = nan(nStarts, 1);
    finalLogLikelihood = nan(nStarts, 1);
    emIterations = nan(nStarts, 1);
    converged = false(nStarts, 1);
    logLikelihoodIsMonotone = false(nStarts, 1);
    hitQFloor = false(nStarts, 1);
    hitQCeiling = false(nStarts, 1);
    hitSigmaBioFloor = false(nStarts, 1);
    hitSigmaBioCeiling = false(nStarts, 1);
    isBestStart = false(nStarts, 1);

    for i = 1:nStarts
        thisStart = emResults.startResults{i};
        initialQF0(i) = thisStart.initialQF0;
        initialQDelta(i) = thisStart.initialQDelta;
        initialSigmaBio(i) = thisStart.initialSigmaBio;
        finalQF0(i) = thisStart.q_f0_em;
        finalQDelta(i) = thisStart.q_delta_em;
        finalSigmaBio(i) = thisStart.sigmaBio_em;
        finalLogLikelihood(i) = thisStart.finalLogLikelihood;
        emIterations(i) = thisStart.emIterations;
        converged(i) = thisStart.converged;
        logLikelihoodIsMonotone(i) = thisStart.logLikelihoodIsMonotone;
        hitQFloor(i) = thisStart.hitQFloor;
        hitQCeiling(i) = thisStart.hitQCeiling;
        hitSigmaBioFloor(i) = thisStart.hitSigmaBioFloor;
        hitSigmaBioCeiling(i) = thisStart.hitSigmaBioCeiling;
    end

    isBestStart(emResults.bestStartIndex) = true;
    summaryTable = table(startMultiplier, initialQF0, initialQDelta, initialSigmaBio, ...
        finalQF0, finalQDelta, finalSigmaBio, finalLogLikelihood, emIterations, ...
        converged, logLikelihoodIsMonotone, hitQFloor, hitQCeiling, ...
        hitSigmaBioFloor, hitSigmaBioCeiling, isBestStart);
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
