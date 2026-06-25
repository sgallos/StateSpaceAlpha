function run_approach1_mfdb_noise_only()
%RUN_APPROACH1_MFDB_NOISE_ONLY Approach 1: MFDB observation noise only.
%
% In this analysis, the SSM observation noise variance is exactly the
% table's mfdb_var column:
%
%   sigma2_obs_k = mfdb_var(k)
%   sigma2_bio   = 0
%
% Two paths differ only in the source/meaning of mfdb_var:
%   Path A: original 120 s MFDB spectral-estimation variance
%   Path B: collapsed 4 x 30 s sub-epoch SEM variance
%
% Each path is fit twice:
%   1) grid-q: q_f0 and q_delta are selected by likelihood grid search
%   2) fixed-q: q_f0 = q_delta = end-to-end-rate heuristic

clearvars -except ans; clc;

%% Step 1: Locations and shared q grid
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outputsDir = fullfile(repoRoot, 'outputs');
figuresDir = fullfile(repoRoot, '4_figures');
addpath(fullfile(repoRoot, '2_state_space_model'));

if ~exist(outputsDir, 'dir')
    mkdir(outputsDir);
end
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

qF0Grid = logspace(-2, 3, 30);
qDeltaGrid = logspace(-2, 3, 30);

paths = struct( ...
    'name', {'PathA_120s', 'PathB_subepoch'}, ...
    'tableFile', { ...
        fullfile(outputsDir, 'alpha_table_for_ssm.csv'), ...
        fullfile(outputsDir, 'alpha_table_for_ssm_subepochs_collapsed.csv')}, ...
    'label', { ...
        'Path A: pure 120s MFDB variance', ...
        'Path B: sub-epoch SEM variance'});

summaryRows = cell(numel(paths) * 2, 1);
summaryIdx = 0;

%% Step 2: Run both paths
for p = 1:numel(paths)
    thisPath = paths(p);
    fprintf('\n============================================================\n');
    fprintf('Approach 1: %s\n', thisPath.label);
    fprintf('Observation noise: sigma2_obs = mfdb_var, sigma2_bio = 0\n');
    fprintf('============================================================\n');

    if ~isfile(thisPath.tableFile)
        warning('Missing table for %s:\n  %s\nSkipping.', thisPath.name, thisPath.tableFile);
        continue;
    end

    alphaTable = readtable(thisPath.tableFile, 'TextType', 'string');

    gridResult = run_one_path_grid_q(alphaTable, qF0Grid, qDeltaGrid, ...
        thisPath, outputsDir, figuresDir);
    summaryIdx = summaryIdx + 1;
    summaryRows{summaryIdx} = make_summary_row(thisPath, "grid-q", gridResult);

    qFixed = compute_q_init(alphaTable.ageYears, alphaTable.alpha_dB);
    fixedResult = run_one_path_fixed_q(alphaTable, qFixed, ...
        thisPath, outputsDir, figuresDir);
    summaryIdx = summaryIdx + 1;
    summaryRows{summaryIdx} = make_summary_row(thisPath, "fixed-q", fixedResult);

    fprintf('\n%s grid-q:  q_f0 = %.4g, q_delta = %.4g, logLik = %.3f, excludes zero = %d\n', ...
        thisPath.name, gridResult.qF0, gridResult.qDelta, gridResult.logLikelihood, ...
        gridResult.excludesZeroAnywhere);
    fprintf('%s fixed-q: q_f0 = %.4g, q_delta = %.4g, logLik = %.3f, excludes zero = %d\n', ...
        thisPath.name, fixedResult.qF0, fixedResult.qDelta, fixedResult.logLikelihood, ...
        fixedResult.excludesZeroAnywhere);
end

summaryRows = summaryRows(1:summaryIdx);
if ~isempty(summaryRows)
    summaryTable = vertcat(summaryRows{:});
    summaryFile = fullfile(outputsDir, 'approach1_mfdb_noise_only_summary.csv');
    writetable(summaryTable, summaryFile);
    fprintf('\nSaved Approach 1 summary to %s\n', summaryFile);
    disp(summaryTable);
else
    warning('No Approach 1 paths were run.');
end

fprintf('\nApproach 1 complete. Figures saved to %s\n', figuresDir);
end

function result = run_one_path_grid_q(alphaTable, qF0Grid, qDeltaGrid, thisPath, outputsDir, figuresDir)
%RUN_ONE_PATH_GRID_Q Grid search over q_f0 and q_delta with sigma2_bio = 0.

nF0 = numel(qF0Grid);
nDelta = numel(qDeltaGrid);
logLikSurface = nan(nF0, nDelta);

for i = 1:nF0
    for j = 1:nDelta
        thisFit = SS_age_diff(alphaTable, ...
            'processNoiseQF0', qF0Grid(i), ...
            'processNoiseQDelta', qDeltaGrid(j), ...
            'biologicalVariance', 0, ...
            'verbose', false);
        logLikSurface(i, j) = thisFit.logLikelihood;
    end
    fprintf('  %s grid row %d/%d complete\n', thisPath.name, i, nF0);
end

[bestLogLik, linIdx] = max(logLikSurface(:));
[iBest, jBest] = ind2sub(size(logLikSurface), linIdx);
bestQF0 = qF0Grid(iBest);
bestQDelta = qDeltaGrid(jBest);

hitGridBoundary = iBest == 1 || iBest == nF0 || jBest == 1 || jBest == nDelta;
if hitGridBoundary
    warning(['%s grid-q maximum is on a q-grid boundary ' ...
        '(q_f0 idx %d/%d, q_delta idx %d/%d).'], ...
        thisPath.name, iBest, nF0, jBest, nDelta);
end

bestFit = SS_age_diff(alphaTable, ...
    'processNoiseQF0', bestQF0, ...
    'processNoiseQDelta', bestQDelta, ...
    'biologicalVariance', 0, ...
    'verbose', false);

traj = sortrows(bestFit.trajectory, 'ageYears');
result = package_fit_result(bestFit, traj, bestQF0, bestQDelta, bestLogLik, ...
    hitGridBoundary, "grid-q");

plot_delta_band(traj, thisPath, "grid-q", bestQF0, bestQDelta, figuresDir);
plot_loglik_heatmap(logLikSurface, qF0Grid, qDeltaGrid, bestQF0, bestQDelta, ...
    thisPath, figuresDir);

save(fullfile(outputsDir, sprintf('approach1_%s_results.mat', thisPath.name)), ...
    'logLikSurface', 'qF0Grid', 'qDeltaGrid', 'bestQF0', 'bestQDelta', ...
    'bestLogLik', 'hitGridBoundary', 'bestFit', 'traj', 'thisPath');
writetable(make_summary_row(thisPath, "grid-q", result), ...
    fullfile(outputsDir, sprintf('approach1_%s_summary.csv', thisPath.name)));
end

function result = run_one_path_fixed_q(alphaTable, qFixed, thisPath, outputsDir, figuresDir)
%RUN_ONE_PATH_FIXED_Q Fit with q fixed by hand and sigma2_bio = 0.

bestFit = SS_age_diff(alphaTable, ...
    'processNoiseQF0', qFixed, ...
    'processNoiseQDelta', qFixed, ...
    'biologicalVariance', 0, ...
    'verbose', false);

traj = sortrows(bestFit.trajectory, 'ageYears');
result = package_fit_result(bestFit, traj, qFixed, qFixed, bestFit.logLikelihood, ...
    false, "fixed-q");

plot_delta_band(traj, thisPath, "fixed-q", qFixed, qFixed, figuresDir);

save(fullfile(outputsDir, sprintf('approach1_%s_fixedq_results.mat', thisPath.name)), ...
    'qFixed', 'bestFit', 'traj', 'thisPath');
writetable(make_summary_row(thisPath, "fixed-q", result), ...
    fullfile(outputsDir, sprintf('approach1_%s_fixedq_summary.csv', thisPath.name)));
end

function result = package_fit_result(fit, traj, qF0, qDelta, logLikelihood, hitGridBoundary, qSelection)
    deltaMean = traj.deltaMean_dB;
    deltaSD = traj.deltaSD_dB;
    deltaLow = deltaMean - 1.96 * deltaSD;
    deltaHigh = deltaMean + 1.96 * deltaSD;
    excludesZeroMask = deltaLow > 0 | deltaHigh < 0;

    [negativeLow, negativeHigh] = get_exclusion_range(traj.ageYears, deltaHigh < 0);
    [positiveLow, positiveHigh] = get_exclusion_range(traj.ageYears, deltaLow > 0);

    result = struct();
    result.qSelection = qSelection;
    result.qF0 = qF0;
    result.qDelta = qDelta;
    result.logLikelihood = logLikelihood;
    result.hitGridBoundary = hitGridBoundary;
    result.meanDelta_dB = mean(deltaMean, 'omitnan');
    result.minDelta_dB = min(deltaMean);
    result.maxDelta_dB = max(deltaMean);
    result.meanDeltaSD_dB = mean(deltaSD, 'omitnan');
    result.maxDeltaSD_dB = max(deltaSD);
    result.excludesZeroAnywhere = any(excludesZeroMask);
    result.negativeExclusionAgeLow = negativeLow;
    result.negativeExclusionAgeHigh = negativeHigh;
    result.positiveExclusionAgeLow = positiveLow;
    result.positiveExclusionAgeHigh = positiveHigh;
    result.fit = fit;
end

function summaryRow = make_summary_row(thisPath, qSelection, result)
    summaryRow = table( ...
        string(thisPath.name), ...
        string(thisPath.label), ...
        string(qSelection), ...
        result.qF0, ...
        result.qDelta, ...
        0, ...
        result.logLikelihood, ...
        result.hitGridBoundary, ...
        result.meanDelta_dB, ...
        result.minDelta_dB, ...
        result.maxDelta_dB, ...
        result.meanDeltaSD_dB, ...
        result.maxDeltaSD_dB, ...
        result.excludesZeroAnywhere, ...
        result.negativeExclusionAgeLow, ...
        result.negativeExclusionAgeHigh, ...
        result.positiveExclusionAgeLow, ...
        result.positiveExclusionAgeHigh, ...
        'VariableNames', {'pathName','pathLabel','qSelection','qF0','qDelta', ...
        'biologicalVariance','logLikelihood','hitGridBoundary','meanDelta_dB', ...
        'minDelta_dB','maxDelta_dB','meanDeltaSD_dB','maxDeltaSD_dB', ...
        'excludesZeroAnywhere','negativeExclusionAgeLow','negativeExclusionAgeHigh', ...
        'positiveExclusionAgeLow','positiveExclusionAgeHigh'});
end

function plot_delta_band(traj, thisPath, qSelection, qF0, qDelta, figuresDir)
    ages = traj.ageYears;
    deltaMean = traj.deltaMean_dB;
    deltaSD = traj.deltaSD_dB;
    deltaLow = deltaMean - 1.96 * deltaSD;
    deltaHigh = deltaMean + 1.96 * deltaSD;

    fig = figure('Color', 'w', 'Position', [100 100 900 520]);
    ax = axes(fig);
    hold(ax, 'on');

    if qSelection == "fixed-q"
        bandColor = [0.70 0.85 0.75];
    else
        bandColor = [0.70 0.80 0.90];
    end

    patch(ax, [ages; flipud(ages)], [deltaLow; flipud(deltaHigh)], ...
        bandColor, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(ax, ages, deltaMean, '-k', 'LineWidth', 2.2);
    yline(ax, 0, '--', 'Color', [0.45 0.45 0.45]);
    xlabel(ax, 'Age (years)');
    ylabel(ax, 'CP - Control \delta(a) (dB)');

    if qSelection == "fixed-q"
        qLine = sprintf('q FIXED = %.3g (imposed smoothness), \\sigma^2_{bio}=0', qF0);
        suffix = "fixedq_delta_band";
    else
        qLine = sprintf('grid-selected q_{f0}=%.3g, q_{\\delta}=%.3g, \\sigma^2_{bio}=0', ...
            qF0, qDelta);
        suffix = "delta_band";
    end

    title(ax, {thisPath.label, ['MFDB observation noise only: ' qLine]}, ...
        'FontWeight', 'normal');
    grid(ax, 'on');
    box(ax, 'off');

    exportgraphics(fig, fullfile(figuresDir, ...
        sprintf('approach1_%s_%s.png', thisPath.name, suffix)), 'Resolution', 200);
    close(fig);
end

function plot_loglik_heatmap(logLikSurface, qF0Grid, qDeltaGrid, bestQF0, bestQDelta, thisPath, figuresDir)
    fig = figure('Color', 'w', 'Position', [100 100 760 640]);
    ax = axes(fig);

    imagesc(ax, log10(qDeltaGrid), log10(qF0Grid), logLikSurface);
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    plot(ax, log10(bestQDelta), log10(bestQF0), 'rp', ...
        'MarkerSize', 16, 'MarkerFaceColor', 'r');
    colorbar(ax);
    xlabel(ax, 'log_{10}(q_\delta)');
    ylabel(ax, 'log_{10}(q_{f0})');
    title(ax, {thisPath.label, ...
        'Log-likelihood over q grid (red star = maximum)'}, ...
        'FontWeight', 'normal');

    exportgraphics(fig, fullfile(figuresDir, ...
        sprintf('approach1_%s_loglik_heatmap.png', thisPath.name)), 'Resolution', 200);
    close(fig);
end

function qInit = compute_q_init(ageYears, alphaDB)
    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange ^ 2) / max(ageSpan ^ 3, eps);
end

function [ageLow, ageHigh] = get_exclusion_range(ageYears, mask)
    if any(mask)
        ageLow = min(ageYears(mask));
        ageHigh = max(ageYears(mask));
    else
        ageLow = NaN;
        ageHigh = NaN;
    end
end
