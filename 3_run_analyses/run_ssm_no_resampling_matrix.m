%% run_ssm_no_resampling_matrix
% No-subject-resampling comparison matrix for the SSM posterior analysis.
%
% Every run reports the same inferential object: the smoother posterior
% trajectory delta(a) = CP - Control alpha power, with its pointwise 95%
% credible band from smoothedCovariance(3,3,k). Subjects are never
% resampled here. The runs differ only in hyperparameter selection.

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

collapsedAlphaTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs_collapsed.csv');
subepochAlphaTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs.csv');
upstreamSubepochAlphaTableFile = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'alpha_table_for_ssm_subepochs.csv');

sigmaBioTrendModel = "linear";
sigmaBioStartForFixedQ = 2.5;

% Profile-Likelihood grid.
qF0Grid = logspace(-2, 3, 30);
qDeltaGrid = logspace(-2, 3, 30);

%% Step 2: Load inputs
if ~isfile(collapsedAlphaTableFile)
    error('Missing collapsed table. Run collapse_subepoch_alpha_table_for_ssm first.');
end
if ~isfile(subepochAlphaTableFile)
    if isfile(upstreamSubepochAlphaTableFile)
        copyfile(upstreamSubepochAlphaTableFile, subepochAlphaTableFile);
        fprintf('Copied sub-epoch table into StateSpaceAlpha outputs:\n  %s\n', subepochAlphaTableFile);
    else
        error('Missing sub-epoch table. Run smoke_test_alpha_extraction_subepochs first.');
    end
end

alphaTable = readtable(collapsedAlphaTableFile, 'TextType', 'string');
subepochTable = readtable(subepochAlphaTableFile, 'TextType', 'string');
[sigmaBioEmpirical, sigmaBioDiagnostics] = estimate_sigmaBio_from_subepochs(subepochTable, ...
    'TrendModel', sigmaBioTrendModel, ...
    'SigmaBioFloor', 0.01, ...
    'Verbose', true);
qInit = compute_q_init(alphaTable.ageYears, alphaTable.alpha_dB);

fprintf('\nNo-resampling matrix setup\n');
fprintf('  qInit: %.6g\n', qInit);
fprintf('  empirical sigmaBio: %.6g dB^2\n', sigmaBioEmpirical);

runResults = struct([]);
trajectoryTables = {};

%% Step 3: Fixed-Heuristic, fixed qInit and fixed empirical sigmaBio
fixedHeuristicName = "Fixed-Heuristic";
fprintf('\n=== %s: fixed qInit, fixed empirical sigmaBio, collapsed 19-row table ===\n', fixedHeuristicName);
fixedHeuristicResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', qInit, ...
    'processNoiseQDelta', qInit, ...
    'biologicalVariance', sigmaBioEmpirical, ...
    'verbose', false);
runResults = make_run_record(fixedHeuristicName, ...
    "fixed qInit", "fixed empirical", "collapsed 19-row", ...
    qInit, qInit, sigmaBioEmpirical, fixedHeuristicResult.logLikelihood, ...
    NaN, true, true, false, fixedHeuristicResult.trajectory);
trajectoryTables{end + 1} = standardize_trajectory(fixedHeuristicName, fixedHeuristicResult.trajectory);

%% Step 4: EM-Smoothness, EM q with fixed empirical sigmaBio
emSmoothnessName = "EM-Smoothness";
fprintf('\n=== %s: EM q, fixed empirical sigmaBio, collapsed 19-row table ===\n', emSmoothnessName);
emSmoothnessResults = SS_age_diff_em(alphaTable, ...
    'fixSigmaBio', true, ...
    'initialSigmaBio', sigmaBioEmpirical, ...
    'qInitMultipliers', [0.01 0.1 1 10 100], ...
    'maxIter', 100, ...
    'tolerance', 1e-5, ...
    'useApproximateMstep', false, ...
    'verbose', false);
emSmoothnessResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', emSmoothnessResults.q_f0_em, ...
    'processNoiseQDelta', emSmoothnessResults.q_delta_em, ...
    'biologicalVariance', sigmaBioEmpirical, ...
    'verbose', false);
runResults(end + 1) = make_run_record(emSmoothnessName, ...
    "EM q", "fixed empirical", "collapsed 19-row", ...
    emSmoothnessResults.q_f0_em, emSmoothnessResults.q_delta_em, ...
    sigmaBioEmpirical, emSmoothnessResult.logLikelihood, ...
    emSmoothnessResults.bestStart.emIterations, emSmoothnessResults.converged, ...
    emSmoothnessResults.bestStart.logLikelihoodIsMonotone, ...
    emSmoothnessResults.hitQFloor || emSmoothnessResults.hitQCeiling, ...
    emSmoothnessResult.trajectory);
trajectoryTables{end + 1} = standardize_trajectory(emSmoothnessName, emSmoothnessResult.trajectory);

%% Step 5: EM-Joint, joint EM on collapsed 19-row table
emJointName = "EM-Joint";
fprintf('\n=== %s: joint EM q + sigmaBio, collapsed 19-row table ===\n', emJointName);
[emJointResults, cSource] = load_or_run_joint_em_collapsed(alphaTable);
emJointResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', emJointResults.q_f0_em, ...
    'processNoiseQDelta', emJointResults.q_delta_em, ...
    'biologicalVariance', emJointResults.sigmaBio_em, ...
    'verbose', false);
runResults(end + 1) = make_run_record(emJointName, ...
    "joint EM q+sigma", "joint EM", "collapsed 19-row", ...
    emJointResults.q_f0_em, emJointResults.q_delta_em, emJointResults.sigmaBio_em, ...
    emJointResult.logLikelihood, ...
    emJointResults.bestStart.emIterations, emJointResults.converged, ...
    emJointResults.bestStart.logLikelihoodIsMonotone, ...
    emJointResults.hitQFloor || emJointResults.hitQCeiling || ...
    emJointResults.hitSigmaBioFloor || emJointResults.hitSigmaBioCeiling, ...
    emJointResult.trajectory);
runResults(end).source = cSource;
trajectoryTables{end + 1} = standardize_trajectory(emJointName, emJointResult.trajectory);

%% Step 6: EM-Joint-Subepoch, joint EM on sub-epoch 76-row table
emJointSubepochName = "EM-Joint-Subepoch";
fprintf('\n=== %s: joint EM q + sigmaBio, sub-epoch 76-row table ===\n', emJointSubepochName);
[emJointSubepochResults, emJointSubepochTrajectory, dSource] = load_or_run_joint_em_subepoch76(subepochTable);
runResults(end + 1) = make_run_record(emJointSubepochName, ...
    "joint EM q+sigma", "joint EM", "sub-epoch 76-row", ...
    emJointSubepochResults.q_f0_em, emJointSubepochResults.q_delta_em, ...
    emJointSubepochResults.sigmaBio_em, emJointSubepochResults.finalLogLikelihood, ...
    emJointSubepochResults.bestStart.emIterations, emJointSubepochResults.converged, ...
    emJointSubepochResults.bestStart.logLikelihoodIsMonotone, ...
    emJointSubepochResults.hitQFloor || emJointSubepochResults.hitQCeiling || ...
    emJointSubepochResults.hitSigmaBioFloor || emJointSubepochResults.hitSigmaBioCeiling, ...
    emJointSubepochTrajectory);
runResults(end).source = dSource;
trajectoryTables{end + 1} = standardize_trajectory(emJointSubepochName, emJointSubepochTrajectory);

%% Step 7: EM-Biological, fixed qInit with EM sigmaBio
emBiologicalName = "EM-Biological";
fprintf('\n=== %s: fixed qInit, EM sigmaBio, collapsed 19-row table ===\n', emBiologicalName);
emBiologicalResults = SS_age_diff_em(alphaTable, ...
    'fixQ', true, ...
    'initialProcessNoiseQF0', qInit, ...
    'initialProcessNoiseQDelta', qInit, ...
    'initialSigmaBio', sigmaBioStartForFixedQ, ...
    'qInitMultipliers', 1, ...
    'maxIter', 100, ...
    'tolerance', 1e-5, ...
    'useApproximateMstep', false, ...
    'verbose', false);
emBiologicalResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', qInit, ...
    'processNoiseQDelta', qInit, ...
    'biologicalVariance', emBiologicalResults.sigmaBio_em, ...
    'verbose', false);
runResults(end + 1) = make_run_record(emBiologicalName, ...
    "fixed qInit", "EM sigma", "collapsed 19-row", ...
    qInit, qInit, emBiologicalResults.sigmaBio_em, emBiologicalResult.logLikelihood, ...
    emBiologicalResults.bestStart.emIterations, emBiologicalResults.converged, ...
    emBiologicalResults.bestStart.logLikelihoodIsMonotone, ...
    emBiologicalResults.hitSigmaBioFloor || emBiologicalResults.hitSigmaBioCeiling, ...
    emBiologicalResult.trajectory);
trajectoryTables{end + 1} = standardize_trajectory(emBiologicalName, emBiologicalResult.trajectory);

%% Step 8: Profile-Likelihood, profile-likelihood q with fixed empirical sigmaBio
profileLikelihoodName = "Profile-Likelihood";
fprintf('\n=== %s: profile-likelihood q, fixed empirical sigmaBio, collapsed 19-row table ===\n', profileLikelihoodName);
[qF0Profile, qDeltaProfile, maxProfileLogLik, logLikSurface] = ...
    profile_likelihood_q(alphaTable, qF0Grid, qDeltaGrid, sigmaBioEmpirical);
profileLikelihoodResult = SS_age_diff(alphaTable, ...
    'processNoiseQF0', qF0Profile, ...
    'processNoiseQDelta', qDeltaProfile, ...
    'biologicalVariance', sigmaBioEmpirical, ...
    'verbose', false);
profileHitBoundary = qF0Profile == qF0Grid(1) || qF0Profile == qF0Grid(end) || ...
    qDeltaProfile == qDeltaGrid(1) || qDeltaProfile == qDeltaGrid(end);
runResults(end + 1) = make_run_record(profileLikelihoodName, ...
    "profile likelihood q", "fixed empirical", "collapsed 19-row", ...
    qF0Profile, qDeltaProfile, sigmaBioEmpirical, profileLikelihoodResult.logLikelihood, ...
    NaN, true, true, profileHitBoundary, profileLikelihoodResult.trajectory);
runResults(end).profileMaxLogLikelihood = maxProfileLogLik;
trajectoryTables{end + 1} = standardize_trajectory(profileLikelihoodName, profileLikelihoodResult.trajectory);

%% Step 9: Save standardized outputs
summaryTable = struct2table(runResults);
allTrajectories = vertcat(trajectoryTables{:});

summaryFile = fullfile(outDir, 'ssm_no_resampling_matrix_summary.csv');
trajectoryFile = fullfile(outDir, 'ssm_no_resampling_matrix_trajectories.csv');
matFile = fullfile(outDir, 'ssm_no_resampling_matrix_results.mat');
profileSurfaceFile = fullfile(outDir, 'ssm_no_resampling_matrix_profile_loglik_surface.csv');

profileSurfaceTable = array2table(logLikSurface, ...
    'VariableNames', matlab.lang.makeUniqueStrings("qDelta_" + string(1:numel(qDeltaGrid))));
profileSurfaceTable.qF0 = qF0Grid(:);
profileSurfaceTable = movevars(profileSurfaceTable, 'qF0', 'Before', 1);

writetable(summaryTable, summaryFile);
writetable(allTrajectories, trajectoryFile);
writetable(profileSurfaceTable, profileSurfaceFile);
save(matFile, 'runResults', 'summaryTable', 'allTrajectories', 'sigmaBioEmpirical', ...
    'sigmaBioDiagnostics', 'qInit', 'qF0Grid', 'qDeltaGrid', 'logLikSurface', ...
    'fixedHeuristicResult', 'emSmoothnessResults', 'emSmoothnessResult', ...
    'emJointResults', 'emJointResult', 'emJointSubepochResults', ...
    'emJointSubepochTrajectory', 'emBiologicalResults', 'emBiologicalResult', ...
    'profileLikelihoodResult');

fprintf('\nSaved no-resampling matrix summary to %s\n', summaryFile);
fprintf('Saved no-resampling matrix trajectories to %s\n', trajectoryFile);
fprintf('Saved profile surface CSV to %s\n', profileSurfaceFile);
fprintf('Saved MAT results to %s\n', matFile);

%% Step 10: Plots
plot_delta_overlay(allTrajectories, fullfile(outDir, 'sensitivity_all_six_overlay.png'));
plot_delta_bands(allTrajectories, fullfile(outDir, 'sensitivity_all_six_bands.png'));
plot_profile_surface(qF0Grid, qDeltaGrid, logLikSurface, qF0Profile, qDeltaProfile, ...
    fullfile(outDir, 'sensitivity_profile_loglik_heatmap.png'));

fprintf('\nNo-resampling matrix complete.\n');
disp(summaryTable(:, {'strategyName','qF0','qDelta','sigmaBio','logLikelihood', ...
    'converged','logLikelihoodIsMonotone','hitBoundary','bandExcludesZero', ...
    'negativeExclusionAgeLow','negativeExclusionAgeHigh'}));

%% Local helpers
function qInit = compute_q_init(ageYears, alphaDB)
    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange ^ 2) / max(ageSpan ^ 3, eps);
end

function [emResults, sourceName] = load_or_run_joint_em_collapsed(alphaTable)
    outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs');
    candidates = [
        string(fullfile(outDir, 'ssm_step4_em_age_difference_results_subepochs_collapsed.mat'))
        string(fullfile(outDir, 'ssm_step4_em_age_difference_results.mat'))
    ];

    emResults = [];
    sourceName = "computed";
    for i = 1:numel(candidates)
        if isfile(candidates(i))
            S = load(candidates(i));
            if isfield(S, 'emResultsCollapsed')
                emResults = S.emResultsCollapsed;
            elseif isfield(S, 'emResults')
                emResults = S.emResults;
            end
            if ~isempty(emResults)
                [~, sourceBase, sourceExt] = fileparts(candidates(i));
                sourceName = "outputs/" + string(sourceBase) + string(sourceExt);
                return;
            end
        end
    end

    emResults = SS_age_diff_em(alphaTable, ...
        'maxIter', 100, ...
        'tolerance', 1e-5, ...
        'qInitMultipliers', [0.01 0.1 1 10 100], ...
        'useApproximateMstep', false, ...
        'verbose', false);
end

function [emResults, subjectTrajectory, sourceName] = load_or_run_joint_em_subepoch76(subepochTable)
    outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'outputs');
    candidate = fullfile(outDir, 'ssm_step4_em_age_difference_results_subepoch76.mat');

    if isfile(candidate)
        S = load(candidate);
        emResults = S.emResults;
        if isfield(S, 'subjectTrajectory')
            subjectTrajectory = S.subjectTrajectory;
        else
            subjectTrajectory = keep_last_subepoch_per_subject(S.rowTrajectory);
        end
        [~, sourceBase, sourceExt] = fileparts(candidate);
        sourceName = "outputs/" + string(sourceBase) + string(sourceExt);
        return;
    end

    emResults = SS_age_diff_em(subepochTable, ...
        'allowRepeatedSubjectRows', true, ...
        'useApproximateMstep', true, ...
        'maxIter', 100, ...
        'tolerance', 1e-5, ...
        'qInitMultipliers', [0.01 0.1 1 10 100], ...
        'verbose', false);
    subjectTrajectory = keep_last_subepoch_per_subject(emResults.trajectory);
    sourceName = "computed";
end

function subjectTrajectory = keep_last_subepoch_per_subject(rowTrajectory)
    [G, subjectIDs] = findgroups(rowTrajectory.subjectID);
    keepMask = false(height(rowTrajectory), 1);
    for i = 1:numel(subjectIDs)
        rows = find(G == i);
        if ismember('subepochIdx', rowTrajectory.Properties.VariableNames)
            [~, orderIdx] = sort(rowTrajectory.subepochIdx(rows));
            rows = rows(orderIdx);
        end
        keepMask(rows(end)) = true;
    end
    subjectTrajectory = sortrows(rowTrajectory(keepMask, :), {'ageYears','subjectID'});
end

function [qF0Best, qDeltaBest, maxLogLik, logLikSurface] = ...
        profile_likelihood_q(alphaTable, qF0Grid, qDeltaGrid, sigmaBioFixed)
    logLikSurface = nan(numel(qF0Grid), numel(qDeltaGrid));
    for i = 1:numel(qF0Grid)
        for j = 1:numel(qDeltaGrid)
            thisFit = SS_age_diff(alphaTable, ...
                'processNoiseQF0', qF0Grid(i), ...
                'processNoiseQDelta', qDeltaGrid(j), ...
                'biologicalVariance', sigmaBioFixed, ...
                'verbose', false);
            logLikSurface(i, j) = thisFit.logLikelihood;
        end
    end

    [maxLogLik, maxIdx] = max(logLikSurface(:));
    [iMax, jMax] = ind2sub(size(logLikSurface), maxIdx);
    qF0Best = qF0Grid(iMax);
    qDeltaBest = qDeltaGrid(jMax);

    fprintf('Profile-likelihood max: q_f0 %.6g | q_delta %.6g | logLik %.6f\n', ...
        qF0Best, qDeltaBest, maxLogLik);
end

function record = make_run_record(strategyName, hyperparameterMethod, sigmaBioMethod, dataStructure, ...
        qF0, qDelta, sigmaBio, logLikelihood, emIterations, converged, ...
        logLikelihoodIsMonotone, hitBoundary, trajectory)
    deltaLow = trajectory.deltaCILow_dB;
    deltaHigh = trajectory.deltaCIHigh_dB;
    ageYears = trajectory.ageYears;
    negativeMask = deltaHigh < 0;
    positiveMask = deltaLow > 0;
    bandExcludesZero = any(negativeMask) || any(positiveMask);

    record = struct();
    record.strategyName = string(strategyName);
    record.hyperparameterMethod = string(hyperparameterMethod);
    record.sigmaBioMethod = string(sigmaBioMethod);
    record.dataStructure = string(dataStructure);
    record.qF0 = qF0;
    record.qDelta = qDelta;
    record.sigmaBio = sigmaBio;
    record.logLikelihood = logLikelihood;
    record.emIterations = emIterations;
    record.converged = logical(converged);
    record.logLikelihoodIsMonotone = logical(logLikelihoodIsMonotone);
    record.hitBoundary = logical(hitBoundary);
    record.meanDelta = mean(trajectory.deltaMean_dB, 'omitnan');
    record.minDelta = min(trajectory.deltaMean_dB);
    record.maxDelta = max(trajectory.deltaMean_dB);
    record.maxBandWidth = max(deltaHigh - deltaLow);
    record.bandExcludesZero = bandExcludesZero;
    record.negativeExclusionAgeLow = min_or_nan(ageYears(negativeMask));
    record.negativeExclusionAgeHigh = max_or_nan(ageYears(negativeMask));
    record.positiveExclusionAgeLow = min_or_nan(ageYears(positiveMask));
    record.positiveExclusionAgeHigh = max_or_nan(ageYears(positiveMask));
    record.source = "";
    record.profileMaxLogLikelihood = NaN;
end

function T = standardize_trajectory(strategyName, trajectory)
    nRows = height(trajectory);
    T = table();
    T.strategyName = repmat(string(strategyName), nRows, 1);
    T.subjectID = string(trajectory.subjectID);
    T.groupLabel = string(trajectory.groupLabel);
    T.ageYears = trajectory.ageYears;
    T.alpha_dB = trajectory.alpha_dB;
    T.observationVariance = trajectory.observationVariance;
    T.deltaMean_dB = trajectory.deltaMean_dB;
    T.deltaSD_dB = trajectory.deltaSD_dB;
    T.deltaCILow_dB = trajectory.deltaCILow_dB;
    T.deltaCIHigh_dB = trajectory.deltaCIHigh_dB;
    T.baselineMean_dB = trajectory.baselineMean_dB;
    T.cpMean_dB = trajectory.cpMean_dB;
end

function val = min_or_nan(x)
    if isempty(x)
        val = NaN;
    else
        val = min(x);
    end
end

function val = max_or_nan(x)
    if isempty(x)
        val = NaN;
    else
        val = max(x);
    end
end

function plot_delta_overlay(allTrajectories, outFile)
    fig = figure('Color', 'w', 'Position', [100 100 1050 560]);
    ax = axes(fig);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    labels = unique(allTrajectories.strategyName, 'stable');
    colors = lines(numel(labels));
    for i = 1:numel(labels)
        rows = allTrajectories.strategyName == labels(i);
        plot(ax, allTrajectories.ageYears(rows), allTrajectories.deltaMean_dB(rows), ...
            '-', 'Color', colors(i, :), 'LineWidth', 2, ...
            'DisplayName', char(labels(i)));
    end
    xlabel(ax, 'Age (years)');
    ylabel(ax, 'CP - Control delta(a) (dB)');
    title(ax, 'Sensitivity matrix: delta(a) posterior means');
    legend(ax, 'Location', 'eastoutside');
    style_axes(ax);
    saveas(fig, outFile);
    fprintf('Saved delta overlay to %s\n', outFile);
end

function plot_delta_bands(allTrajectories, outFile)
    labels = unique(allTrajectories.strategyName, 'stable');
    fig = figure('Color', 'w', 'Position', [100 100 1200 780]);
    tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:numel(labels)
        ax = nexttile;
        rows = allTrajectories.strategyName == labels(i);
        ageYears = allTrajectories.ageYears(rows);
        deltaMean = allTrajectories.deltaMean_dB(rows);
        ciLow = allTrajectories.deltaCILow_dB(rows);
        ciHigh = allTrajectories.deltaCIHigh_dB(rows);
        hold(ax, 'on');
        patch(ax, [ageYears; flipud(ageYears)], [ciLow; flipud(ciHigh)], ...
            [0.72 0.82 0.92], 'EdgeColor', 'none', 'FaceAlpha', 0.55);
        plot(ax, ageYears, deltaMean, '-k', 'LineWidth', 2);
        yline(ax, 0, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
        title(ax, labels(i));
        xlabel(ax, 'Age (years)');
        ylabel(ax, 'delta(a) dB');
        style_axes(ax);
    end
    saveas(fig, outFile);
    fprintf('Saved delta band panels to %s\n', outFile);
end

function plot_profile_surface(qF0Grid, qDeltaGrid, logLikSurface, qF0Best, qDeltaBest, outFile)
    fig = figure('Color', 'w', 'Position', [100 100 790 650]);
    ax = axes(fig);
    imagesc(ax, log10(qDeltaGrid), log10(qF0Grid), logLikSurface);
    axis(ax, 'xy');
    hold(ax, 'on');
    plot(ax, log10(qDeltaBest), log10(qF0Best), 'wo', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.4);
    xlabel(ax, 'log10 q_delta');
    ylabel(ax, 'log10 q_f0');
    title(ax, 'Profile-Likelihood smoothness surface');
    cb = colorbar(ax);
    ylabel(cb, 'Log-likelihood');
    style_axes(ax);
    saveas(fig, outFile);
    fprintf('Saved profile likelihood heatmap to %s\n', outFile);
end

function style_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
