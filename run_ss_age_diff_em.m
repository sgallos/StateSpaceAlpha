%% run_ss_age_diff_em
% Step 4 entry point for estimating the IWP smoothness hyperparameters by EM.
%
% This script loads the 19-row alpha table from the smoke test, estimates
% q_f0, q_delta, and biological variance by multi-start EM, and writes the EM-selected
% baseline/control, CP, and CP-Control difference trajectories.

clear; clc;

%% Step 1: Locate the alpha table
repoRoot = fileparts(mfilename('fullpath'));
outDir = fullfile(repoRoot, 'outputs');

% Choose which 19-row alpha table to fit.
%   alpha_table_for_ssm.csv = original 120 s MFDB table
%   alpha_table_for_ssm_subepochs_collapsed.csv = collapsed 4 x 30 s table
alphaTableFileName = 'alpha_table_for_ssm.csv';
alphaTableFile = fullfile(outDir, alphaTableFileName);

if ~isfile(alphaTableFile)
    error(['Alpha table not found:\n%s\n\n' ...
        'Run this first:\n' ...
        '  cd(''%s'')\n' ...
        '  smoke_test_alpha_extraction'], alphaTableFile, repoRoot);
end

alphaTable = readtable(alphaTableFile, 'TextType', 'string');
fprintf('Loaded alpha table: %s\n', alphaTableFile);
disp(alphaTable(:, {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var'}));

%% Step 2: Run multi-start EM
emOptions = struct();
emOptions.maxIter = 100;
emOptions.tolerance = 1e-5;
emOptions.qInitMultipliers = [0.01 0.1 1 10 100];
emOptions.useApproximateMstep = false;
emOptions.initialSigmaBio = 2.451431;  % from run_ssm_residual_sanity_check

emResults = SS_age_diff_em(alphaTable, ...
    'maxIter', emOptions.maxIter, ...
    'tolerance', emOptions.tolerance, ...
    'qInitMultipliers', emOptions.qInitMultipliers, ...
    'initialSigmaBio', emOptions.initialSigmaBio, ...
    'useApproximateMstep', emOptions.useApproximateMstep, ...
    'verbose', true);

bestResult = emResults.bestStart;
bestTrajectory = emResults.trajectory;

if emResults.hitQFloor || emResults.hitQCeiling || emResults.hitSigmaBioFloor || emResults.hitSigmaBioCeiling
    warning(['The best EM solution hit a q boundary. This usually means the current ' ...
        'observation-variance model is driving the likelihood toward an extreme ' ...
        'smoothness setting. Inspect the Step 4 diagnostics before interpreting the trajectory.']);
end

%% Step 3: Save Step 4 outputs
step4MatFile = fullfile(outDir, 'ssm_step4_em_age_difference_results.mat');
step4CsvFile = fullfile(outDir, 'ssm_step4_em_age_difference_trajectory.csv');
step4SummaryFile = fullfile(outDir, 'ssm_step4_em_multistart_summary.csv');

multiStartSummary = build_multistart_summary(emResults);

writetable(bestTrajectory, step4CsvFile);
writetable(multiStartSummary, step4SummaryFile);
save(step4MatFile, 'alphaTable', 'emResults', 'bestResult', ...
    'bestTrajectory', 'multiStartSummary', 'emOptions', ...
    'alphaTableFileName', 'alphaTableFile');

fprintf('\nSaved Step 4 trajectory CSV to %s\n', step4CsvFile);
fprintf('Saved Step 4 multi-start summary to %s\n', step4SummaryFile);
fprintf('Saved Step 4 MAT results to %s\n', step4MatFile);

%% Step 4: Figure A, EM-selected CP-Control difference
figDelta = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figDelta);
hold(ax, 'on');

hBand = plot_ci_band(ax, bestTrajectory.ageYears, ...
    bestTrajectory.deltaCILow_dB, bestTrajectory.deltaCIHigh_dB, ...
    [0.78 0.80 0.86], 0.42);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hDelta = plot(ax, bestTrajectory.ageYears, bestTrajectory.deltaMean_dB, ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, sprintf('Step 4 EM-selected CP-Control alpha-power difference (qF0 %.3g, qDelta %.3g, sigmaBio %.3g)', ...
    emResults.q_f0_em, emResults.q_delta_em, emResults.sigmaBio_em));
legend(ax, [hBand, hDelta, hZero], {'Delta 95% CI', 'Delta mean', 'Zero line'}, ...
    'Location', 'best');
style_ssm_axes(ax);

deltaPlotFile = fullfile(outDir, 'ssm_step4_em_delta_primary.png');
saveas(figDelta, deltaPlotFile);
fprintf('Saved Step 4 EM delta plot to %s\n', deltaPlotFile);

%% Step 5: Figure B, EM log-likelihood histories
figLogLik = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figLogLik);
hold(ax, 'on');

for i = 1:numel(emResults.startResults)
    thisStart = emResults.startResults{i};
    plot(ax, 1:numel(thisStart.logLikHistory), thisStart.logLikHistory, ...
        '-o', 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'DisplayName', sprintf('start %.2g', emResults.qInitMultipliers(i)));
end

xlabel(ax, 'EM iteration');
ylabel(ax, 'Log-likelihood');
title(ax, 'Step 4 EM convergence diagnostic');
legend(ax, 'Location', 'best');
style_ssm_axes(ax);

logLikPlotFile = fullfile(outDir, 'ssm_step4_em_loglik_history.png');
saveas(figLogLik, logLikPlotFile);
fprintf('Saved Step 4 EM log-likelihood plot to %s\n', logLikPlotFile);

%% Step 6: Figure C, hyperparameter histories
figQ = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figQ);
hold(ax, 'on');

for i = 1:numel(emResults.startResults)
    thisStart = emResults.startResults{i};
    iterVals = 1:size(thisStart.qHistory, 1);
    semilogy(ax, iterVals, thisStart.qHistory(:, 1), '-', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q f0 start %.2g', emResults.qInitMultipliers(i)));
    semilogy(ax, iterVals, thisStart.qHistory(:, 2), '--', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q delta start %.2g', emResults.qInitMultipliers(i)));
    semilogy(ax, iterVals, thisStart.sigmaBioHistory, ':', ...
        'LineWidth', 2.0, 'DisplayName', sprintf('sigma bio start %.2g', emResults.qInitMultipliers(i)));
end

xlabel(ax, 'EM iteration');
ylabel(ax, 'Hyperparameter value');
title(ax, 'Step 4 EM hyperparameter histories');
legend(ax, 'Location', 'eastoutside');
set(ax, 'YScale', 'log');
style_ssm_axes(ax);

qHistoryPlotFile = fullfile(outDir, 'ssm_step4_em_hyperparameter_history.png');
saveas(figQ, qHistoryPlotFile);
fprintf('Saved Step 4 EM hyperparameter-history plot to %s\n', qHistoryPlotFile);

%% Step 7: Optional comparison against Step 3 fixed-q sensitivity
step3MatFile = fullfile(outDir, 'ssm_step3_age_difference_results.mat');
if isfile(step3MatFile)
    S = load(step3MatFile, 'qScalesF0', 'qScalesDelta', 'sensitivityResults');

    figCompare = figure('Color', 'w', 'Position', [100 100 1120 590]);
    ax = axes(figCompare);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);

    for i = 1:numel(S.qScalesF0)
        for j = 1:numel(S.qScalesDelta)
            thisResult = S.sensitivityResults{i, j};
            plot(ax, thisResult.trajectory.ageYears, ...
                thisResult.trajectory.deltaMean_dB, '-', ...
                'Color', [0.70 0.70 0.70], 'LineWidth', 1.1, ...
                'HandleVisibility', 'off');
        end
    end

    hEM = plot(ax, bestTrajectory.ageYears, bestTrajectory.deltaMean_dB, ...
        '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 3.0);
    xlabel(ax, 'Age (years)');
    ylabel(ax, 'CP - Control alpha power difference (dB)');
    title(ax, 'Step 4 EM trajectory over Step 3 fixed-q sensitivity curves');
    legend(ax, hEM, {'EM-selected delta'}, 'Location', 'best');
    style_ssm_axes(ax);

    comparePlotFile = fullfile(outDir, 'ssm_step4_em_vs_step3_sensitivity.png');
    saveas(figCompare, comparePlotFile);
    fprintf('Saved Step 4 EM vs Step 3 comparison plot to %s\n', comparePlotFile);
end

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
