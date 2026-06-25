function run_power_ssm_scalar_r()
%RUN_POWER_SSM_SCALAR_R Power-band SSM with one EM-fit observation variance.
%
% This workflow estimates Control, CP, and CP-Control group-difference
% power trajectories using a scalar observation-noise variance r:
%
%   R_k = r
%
% MFDB variance and sigmaBio are not used in this workflow.

%% ========================================================================
%% USER SETTINGS -- edit only this block
%% ========================================================================

% --- Which band to run ---
bandName = "alpha";        % "delta" | "theta" | "alpha"
bandLowHz = 8;
bandHighHz = 12;

% --- Input data ---
inputTable = "alpha_table_for_ssm_subepochs_collapsed.csv";
powerColumn = "alpha_dB";

% --- Observation noise mode ---
observationVarianceMode = "scalar_r";  % this workflow expects "scalar_r"

% --- q/state-noise treatment ---
qMode = "all";              % "all" | "EM-Joint" | "FixedHeuristicQ" | "FixedSmoothQ" | "Profile-q"
qFixedValue = [];           % used for FixedSmoothQ; [] means 0.1 * qInit
qF0Grid = logspace(-2, 3, 30);
qGroupDifferenceGrid = logspace(-2, 3, 30);

% --- EM settings ---
emMaxIterations = 100;
emTolerance = 1e-5;
qInitMultipliers = [0.01 0.1 1 10 100];
rFloor = 1e-4;
rCeiling = 1e4;

% --- Output ---
outputPrefix = "power_ssm";

%% ========================================================================
%% END USER SETTINGS -- do not edit below
%% ========================================================================

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

bandConfig = struct();
bandConfig.bandName = bandName;
bandConfig.bandLowHz = bandLowHz;
bandConfig.bandHighHz = bandHighHz;
bandConfig.inputTable = inputTable;
bandConfig.powerColumn = powerColumn;
bandConfig.observationVarianceMode = observationVarianceMode;
bandConfig.outputPrefix = outputPrefix;

inputTablePath = resolve_input_table_path(inputTable, outputsDir);
rawTable = readtable(inputTablePath, 'TextType', 'string');
ssmTable = build_ssm_power_table(rawTable, powerColumn, observationVarianceMode);
qInit = compute_q_init(ssmTable.ageYears, ssmTable.alpha_dB);

if isempty(qFixedValue)
    qFixedValue = 0.1 * qInit;
end
bandConfig.qInit = qInit;
bandConfig.qFixedValue = qFixedValue;

fprintf('\nPower SSM scalar-r workflow\n');
fprintf('  band: %s [%.3g, %.3g] Hz\n', bandName, bandLowHz, bandHighHz);
fprintf('  input: %s\n', inputTablePath);
fprintf('  power column: %s\n', powerColumn);
fprintf('  qInit: %.6g | FixedSmoothQ q: %.6g\n', qInit, qFixedValue);
fprintf('  observation variance mode: %s\n', observationVarianceMode);

versionResults = struct([]);

if should_run(qMode, "EM-Joint")
    result = run_em_joint_version(ssmTable, bandConfig, ...
        emMaxIterations, emTolerance, qInitMultipliers, rFloor, rCeiling);
    versionResults = append_version_result(versionResults, result);
end

if should_run(qMode, "FixedHeuristicQ")
    result = run_fixed_q_version(ssmTable, bandConfig, ...
        "FixedHeuristicQ", "Fixed heuristic q", qInit, ...
        emMaxIterations, emTolerance, rFloor, rCeiling);
    versionResults = append_version_result(versionResults, result);
end

if should_run(qMode, "FixedSmoothQ")
    result = run_fixed_q_version(ssmTable, bandConfig, ...
        "FixedSmoothQ", "Fixed smooth q", qFixedValue, ...
        emMaxIterations, emTolerance, rFloor, rCeiling);
    versionResults = append_version_result(versionResults, result);
end

if should_run(qMode, "Profile-q")
    profileR = choose_profile_r(versionResults);
    result = run_profile_q_version(ssmTable, bandConfig, profileR, ...
        qF0Grid, qGroupDifferenceGrid, rFloor, rCeiling);
    versionResults = append_version_result(versionResults, result);
end

if isempty(versionResults)
    error('No q versions were run. Check qMode.');
end

[summaryTable, trajectoryTable] = build_output_tables(versionResults, bandConfig);
displayVersionName = select_display_version(summaryTable);

outputBase = sprintf('%s_%s', outputPrefix, lower(bandName));
allVersionsFile = fullfile(outputsDir, sprintf('%s_all_versions.mat', outputBase));
summaryFile = fullfile(outputsDir, sprintf('%s_summary.csv', outputBase));
trajectoryFile = fullfile(outputsDir, sprintf('%s_trajectories.csv', outputBase));

writetable(summaryTable, summaryFile);
writetable(trajectoryTable, trajectoryFile);
save(allVersionsFile, 'versionResults', 'summaryTable', 'trajectoryTable', ...
    'bandConfig', 'displayVersionName');

plot_three_trajectories(versionResults, displayVersionName, bandConfig, figuresDir);
plot_group_difference(versionResults, displayVersionName, bandConfig, figuresDir);
plot_q_version_comparison(versionResults, bandConfig, figuresDir);
plot_em_diagnostics(versionResults, bandConfig, figuresDir);
plot_profile_q_heatmap(versionResults, bandConfig, figuresDir);

fprintf('\nSaved scalar-r power SSM outputs:\n');
fprintf('  %s\n', allVersionsFile);
fprintf('  %s\n', summaryFile);
fprintf('  %s\n', trajectoryFile);
fprintf('Display version for single-version figures: %s\n', displayVersionName);
disp(summaryTable(:, {'versionName','qF0','qGroupDifference','r','logLikelihood', ...
    'converged','logLikelihoodIsMonotone','hitBoundary','groupDifferenceExcludesZero'}));
end

function inputTablePath = resolve_input_table_path(inputTable, outputsDir)
    inputTable = string(inputTable);
    if isfile(inputTable)
        inputTablePath = char(inputTable);
    else
        inputTablePath = fullfile(outputsDir, char(inputTable));
    end

    if ~isfile(inputTablePath)
        error('Input table not found: %s', inputTablePath);
    end
end

function ssmTable = build_ssm_power_table(rawTable, powerColumn, observationVarianceMode)
    powerColumn = string(powerColumn);
    requiredColumns = ["subjectID", "groupLabel", "ageYears", powerColumn];
    missingColumns = setdiff(requiredColumns, string(rawTable.Properties.VariableNames));
    if ~isempty(missingColumns)
        error('Input table is missing required columns: %s', strjoin(missingColumns, ', '));
    end

    ssmTable = table();
    ssmTable.subjectID = string(rawTable.subjectID);
    ssmTable.groupLabel = string(rawTable.groupLabel);
    ssmTable.ageYears = double(rawTable.ageYears);
    ssmTable.alpha_dB = double(rawTable.(powerColumn));
    ssmTable.mfdb_var = zeros(height(rawTable), 1);

    if observationVarianceMode == "mfdb_plus_bio" && all(ssmTable.mfdb_var == 0)
        error('mfdb_var is all zeros but mode is mfdb_plus_bio. Use scalar_r mode or provide real mfdb_var.');
    end
end

function tf = should_run(qMode, versionName)
    qMode = string(qMode);
    versionName = string(versionName);
    tf = qMode == "all" || qMode == versionName;
end

function versionResults = append_version_result(versionResults, result)
    if isempty(versionResults)
        versionResults = result;
    else
        versionResults(end + 1) = result;
    end
end

function result = run_em_joint_version(ssmTable, bandConfig, maxIter, tolerance, qInitMultipliers, rFloor, rCeiling)
    versionName = "EM-Joint";
    fprintf('\n=== %s ===\n', versionName);
    emResult = SS_age_diff_em(ssmTable, ...
        'observationVarianceMode', 'scalar_r', ...
        'maxIter', maxIter, ...
        'tolerance', tolerance, ...
        'qInitMultipliers', qInitMultipliers, ...
        'rFloor', rFloor, ...
        'rCeiling', rCeiling, ...
        'verbose', false);

    fit = SS_age_diff(ssmTable, ...
        'observationVarianceMode', 'scalar_r', ...
        'observationVarianceR', emResult.r_em, ...
        'processNoiseQF0', emResult.q_f0_em, ...
        'processNoiseQDelta', emResult.q_delta_em, ...
        'verbose', false);

    result = package_version_result(versionName, "EM fit q and r", bandConfig, ...
        emResult.q_f0_em, emResult.q_delta_em, emResult.r_em, fit, emResult, ...
        emResult.hitQFloor || emResult.hitQCeiling || emResult.hitRFloor || emResult.hitRCeiling);
end

function result = run_fixed_q_version(ssmTable, bandConfig, versionName, versionLabel, qValue, maxIter, tolerance, rFloor, rCeiling)
    fprintf('\n=== %s ===\n', versionName);
    emResult = SS_age_diff_em(ssmTable, ...
        'observationVarianceMode', 'scalar_r', ...
        'fixQ', true, ...
        'initialProcessNoiseQF0', qValue, ...
        'initialProcessNoiseQDelta', qValue, ...
        'qInitMultipliers', 1, ...
        'maxIter', maxIter, ...
        'tolerance', tolerance, ...
        'rFloor', rFloor, ...
        'rCeiling', rCeiling, ...
        'verbose', false);

    fit = SS_age_diff(ssmTable, ...
        'observationVarianceMode', 'scalar_r', ...
        'observationVarianceR', emResult.r_em, ...
        'processNoiseQF0', qValue, ...
        'processNoiseQDelta', qValue, ...
        'verbose', false);

    result = package_version_result(versionName, versionLabel, bandConfig, ...
        qValue, qValue, emResult.r_em, fit, emResult, ...
        emResult.hitRFloor || emResult.hitRCeiling);
end

function profileR = choose_profile_r(versionResults)
    profileR = NaN;
    for i = 1:numel(versionResults)
        if versionResults(i).versionName == "EM-Joint" && ...
                versionResults(i).converged && versionResults(i).logLikelihoodIsMonotone && ...
                isfinite(versionResults(i).r)
            profileR = versionResults(i).r;
            fprintf('Profile-q using r from clean EM-Joint: %.6g\n', profileR);
            return;
        end
    end

    for i = 1:numel(versionResults)
        if versionResults(i).versionName == "FixedHeuristicQ" && isfinite(versionResults(i).r)
            profileR = versionResults(i).r;
            fprintf('Profile-q using r from FixedHeuristicQ fallback: %.6g\n', profileR);
            return;
        end
    end

    for i = 1:numel(versionResults)
        if isfinite(versionResults(i).r)
            profileR = versionResults(i).r;
            fprintf('Profile-q using r from %s fallback: %.6g\n', ...
                versionResults(i).versionName, profileR);
            return;
        end
    end

    error('No finite r is available for Profile-q.');
end

function result = run_profile_q_version(ssmTable, bandConfig, profileR, qF0Grid, qGroupDifferenceGrid, rFloor, rCeiling)
    versionName = "Profile-q";
    fprintf('\n=== %s ===\n', versionName);
    profileR = min(max(profileR, rFloor), rCeiling);
    nF0 = numel(qF0Grid);
    nGroupDifference = numel(qGroupDifferenceGrid);
    logLikSurface = nan(nF0, nGroupDifference);

    for i = 1:nF0
        for j = 1:nGroupDifference
            fit = SS_age_diff(ssmTable, ...
                'observationVarianceMode', 'scalar_r', ...
                'observationVarianceR', profileR, ...
                'processNoiseQF0', qF0Grid(i), ...
                'processNoiseQDelta', qGroupDifferenceGrid(j), ...
                'verbose', false);
            logLikSurface(i, j) = fit.logLikelihood;
        end
        fprintf('  Profile-q row %d/%d complete\n', i, nF0);
    end

    [bestLogLikelihood, linearIndex] = max(logLikSurface(:));
    [bestF0Index, bestGroupDifferenceIndex] = ind2sub(size(logLikSurface), linearIndex);
    bestQF0 = qF0Grid(bestF0Index);
    bestQGroupDifference = qGroupDifferenceGrid(bestGroupDifferenceIndex);
    hitBoundary = bestF0Index == 1 || bestF0Index == nF0 || ...
        bestGroupDifferenceIndex == 1 || bestGroupDifferenceIndex == nGroupDifference;

    fit = SS_age_diff(ssmTable, ...
        'observationVarianceMode', 'scalar_r', ...
        'observationVarianceR', profileR, ...
        'processNoiseQF0', bestQF0, ...
        'processNoiseQDelta', bestQGroupDifference, ...
        'verbose', false);

    result = package_version_result(versionName, "Profile likelihood q, fixed r", bandConfig, ...
        bestQF0, bestQGroupDifference, profileR, fit, [], hitBoundary);
    result.profileR = profileR;
    result.profileLogLikSurface = logLikSurface;
    result.profileQF0Grid = qF0Grid;
    result.profileQGroupDifferenceGrid = qGroupDifferenceGrid;
    result.profileMaxLogLikelihood = bestLogLikelihood;
end

function result = package_version_result(versionName, versionLabel, bandConfig, qF0, qGroupDifference, r, fit, emResult, hitBoundary)
    traj = sortrows(fit.trajectory, 'ageYears');
    groupDifferenceLow = traj.deltaCILow_dB;
    groupDifferenceHigh = traj.deltaCIHigh_dB;
    negativeMask = groupDifferenceHigh < 0;
    positiveMask = groupDifferenceLow > 0;

    result = struct();
    result.versionName = string(versionName);
    result.versionLabel = string(versionLabel);
    result.bandName = string(bandConfig.bandName);
    result.qF0 = qF0;
    result.qGroupDifference = qGroupDifference;
    result.r = r;
    result.logLikelihood = fit.logLikelihood;
    result.fit = fit;
    result.emResult = emResult;
    result.trajectory = traj;
    result.converged = false;
    result.anyConverged = false;
    result.logLikelihoodIsMonotone = true;
    result.emIterations = NaN;
    result.hitBoundary = hitBoundary;
    result.profileR = NaN;
    result.profileLogLikSurface = [];
    result.profileQF0Grid = [];
    result.profileQGroupDifferenceGrid = [];
    result.profileMaxLogLikelihood = NaN;

    if ~isempty(emResult)
        result.converged = emResult.converged;
        result.anyConverged = emResult.anyConverged;
        result.logLikelihoodIsMonotone = emResult.bestStart.logLikelihoodIsMonotone;
        result.emIterations = emResult.bestStart.emIterations;
    end

    result.groupDifferenceMean_dB = mean(traj.deltaMean_dB, 'omitnan');
    result.groupDifferenceMin_dB = min(traj.deltaMean_dB);
    result.groupDifferenceMax_dB = max(traj.deltaMean_dB);
    result.groupDifferenceMeanSD_dB = mean(traj.deltaSD_dB, 'omitnan');
    result.groupDifferenceExcludesZero = any(negativeMask | positiveMask);
    [result.negativeExclusionAgeLow, result.negativeExclusionAgeHigh] = ...
        get_age_range(traj.ageYears, negativeMask);
    [result.positiveExclusionAgeLow, result.positiveExclusionAgeHigh] = ...
        get_age_range(traj.ageYears, positiveMask);

    fprintf('%s: q_f0 %.6g | q_groupDifference %.6g | r %.6g | logLik %.6f\n', ...
        result.versionName, result.qF0, result.qGroupDifference, result.r, result.logLikelihood);
end

function [summaryTable, trajectoryTable] = build_output_tables(versionResults, bandConfig)
    summaryRows = cell(numel(versionResults), 1);
    trajectoryRows = cell(numel(versionResults), 1);

    for i = 1:numel(versionResults)
        r = versionResults(i);
        summaryRows{i} = table( ...
            string(r.bandName), ...
            string(r.versionName), ...
            string(r.versionLabel), ...
            r.qF0, ...
            r.qGroupDifference, ...
            r.r, ...
            r.logLikelihood, ...
            r.converged, ...
            r.anyConverged, ...
            r.logLikelihoodIsMonotone, ...
            r.emIterations, ...
            r.hitBoundary, ...
            r.groupDifferenceMean_dB, ...
            r.groupDifferenceMin_dB, ...
            r.groupDifferenceMax_dB, ...
            r.groupDifferenceMeanSD_dB, ...
            r.groupDifferenceExcludesZero, ...
            r.negativeExclusionAgeLow, ...
            r.negativeExclusionAgeHigh, ...
            r.positiveExclusionAgeLow, ...
            r.positiveExclusionAgeHigh, ...
            'VariableNames', {'bandName','versionName','versionLabel','qF0', ...
            'qGroupDifference','r','logLikelihood','converged','anyConverged', ...
            'logLikelihoodIsMonotone','emIterations','hitBoundary', ...
            'groupDifferenceMean_dB','groupDifferenceMin_dB','groupDifferenceMax_dB', ...
            'groupDifferenceMeanSD_dB','groupDifferenceExcludesZero', ...
            'negativeExclusionAgeLow','negativeExclusionAgeHigh', ...
            'positiveExclusionAgeLow','positiveExclusionAgeHigh'});

        thisTrajectory = r.trajectory;
        thisTrajectory.bandName = repmat(string(bandConfig.bandName), height(thisTrajectory), 1);
        thisTrajectory.versionName = repmat(string(r.versionName), height(thisTrajectory), 1);
        thisTrajectory = movevars(thisTrajectory, {'bandName','versionName'}, 'Before', 1);
        trajectoryRows{i} = thisTrajectory;
    end

    summaryTable = vertcat(summaryRows{:});
    trajectoryTable = vertcat(trajectoryRows{:});
end

function displayVersionName = select_display_version(summaryTable)
    cleanEM = summaryTable.versionName == "EM-Joint" & summaryTable.converged & ...
        summaryTable.logLikelihoodIsMonotone & ~summaryTable.hitBoundary;
    if any(cleanEM)
        displayVersionName = "EM-Joint";
        return;
    end

    if any(summaryTable.versionName == "FixedHeuristicQ")
        displayVersionName = "FixedHeuristicQ";
    else
        displayVersionName = summaryTable.versionName(1);
    end
end

function plot_three_trajectories(versionResults, displayVersionName, bandConfig, figuresDir)
    result = get_version_result(versionResults, displayVersionName);
    traj = result.trajectory;
    ages = traj.ageYears;
    isControl = traj.groupIndicator == 0;
    isCP = traj.groupIndicator == 1;

    fig = figure('Color', 'w', 'Position', [80 80 980 900]);
    tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl);
    hold(ax, 'on');
    plot_band(ax, ages, traj.baselineMean_dB, traj.baselineCILow_dB, traj.baselineCIHigh_dB, [0.55 0.68 0.92]);
    scatter(ax, traj.ageYears(isControl), traj.alpha_dB(isControl), 58, [0.12 0.32 0.78], ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    ylabel(ax, sprintf('%s power (dB)', title_case(bandConfig.bandName)));
    title(ax, sprintf('Control trajectory f_0(a) -- %s', result.versionName), 'FontWeight', 'normal');
    style_axis(ax);

    ax = nexttile(tl);
    hold(ax, 'on');
    plot_band(ax, ages, traj.cpMean_dB, traj.cpCILow_dB, traj.cpCIHigh_dB, [0.92 0.62 0.62]);
    scatter(ax, traj.ageYears(isCP), traj.alpha_dB(isCP), 58, [0.78 0.16 0.16], ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    ylabel(ax, sprintf('%s power (dB)', title_case(bandConfig.bandName)));
    title(ax, 'CP trajectory f_0(a) + group difference', 'FontWeight', 'normal');
    style_axis(ax);

    ax = nexttile(tl);
    hold(ax, 'on');
    plot_band(ax, ages, traj.deltaMean_dB, traj.deltaCILow_dB, traj.deltaCIHigh_dB, [0.74 0.74 0.74]);
    cpDifferenceDots = traj.alpha_dB(isCP) - traj.baselineMean_dB(isCP);
    scatter(ax, traj.ageYears(isCP), cpDifferenceDots, 52, [0.45 0.45 0.45], ...
        'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
    yline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);
    xlabel(ax, 'Age (years)');
    ylabel(ax, 'CP - Control (dB)');
    title(ax, 'CP-Control group difference; dots are CP observations minus fitted Control mean', ...
        'FontWeight', 'normal');
    style_axis(ax);

    title(tl, sprintf('%s-band power SSM trajectories (%s)', ...
        title_case(bandConfig.bandName), result.versionName), ...
        'FontWeight', 'normal');

    outFile = fullfile(figuresDir, sprintf('%s_%s_three_trajectories.png', ...
        bandConfig.outputPrefix, lower(bandConfig.bandName)));
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
end

function plot_group_difference(versionResults, displayVersionName, bandConfig, figuresDir)
    result = get_version_result(versionResults, displayVersionName);
    traj = result.trajectory;

    fig = figure('Color', 'w', 'Position', [100 100 900 540]);
    ax = axes(fig);
    hold(ax, 'on');
    plot_band(ax, traj.ageYears, traj.deltaMean_dB, traj.deltaCILow_dB, traj.deltaCIHigh_dB, [0.70 0.80 0.90]);
    yline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);
    xlabel(ax, 'Age (years)');
    ylabel(ax, sprintf('CP - Control %s power (dB)', title_case(bandConfig.bandName)));
    title(ax, {sprintf('%s-band CP-Control group difference', title_case(bandConfig.bandName)), ...
        sprintf('%s: q_{f0}=%.3g, q_{group}=%.3g, r=%.3g dB^2', ...
        result.versionName, result.qF0, result.qGroupDifference, result.r)}, ...
        'FontWeight', 'normal');
    style_axis(ax);

    outFile = fullfile(figuresDir, sprintf('%s_%s_group_difference.png', ...
        bandConfig.outputPrefix, lower(bandConfig.bandName)));
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
end

function plot_q_version_comparison(versionResults, bandConfig, figuresDir)
    fig = figure('Color', 'w', 'Position', [60 60 1300 820]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    colors = lines(max(numel(versionResults), 4));

    for i = 1:numel(versionResults)
        ax = nexttile(tl);
        hold(ax, 'on');
        traj = versionResults(i).trajectory;
        plot_band(ax, traj.ageYears, traj.deltaMean_dB, traj.deltaCILow_dB, ...
            traj.deltaCIHigh_dB, colors(i, :));
        yline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);
        xlabel(ax, 'Age (years)');
        ylabel(ax, 'CP - Control (dB)');
        title(ax, sprintf('%s | r=%.3g, q=(%.3g, %.3g)', ...
            versionResults(i).versionName, versionResults(i).r, ...
            versionResults(i).qF0, versionResults(i).qGroupDifference), ...
            'FontWeight', 'normal');
        style_axis(ax);
    end

    for i = (numel(versionResults) + 1):4
        ax = nexttile(tl);
        axis(ax, 'off');
    end

    title(tl, sprintf('%s-band q-version comparison, scalar-r observation noise', ...
        title_case(bandConfig.bandName)), 'FontWeight', 'normal');
    outFile = fullfile(figuresDir, sprintf('%s_%s_q_versions.png', ...
        bandConfig.outputPrefix, lower(bandConfig.bandName)));
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
end

function plot_em_diagnostics(versionResults, bandConfig, figuresDir)
    emMask = arrayfun(@(r) ~isempty(r.emResult), versionResults);
    emResults = versionResults(emMask);
    if isempty(emResults)
        return;
    end

    fig = figure('Color', 'w', 'Position', [60 60 420 * numel(emResults) 780]);
    tl = tiledlayout(fig, 3, numel(emResults), 'TileSpacing', 'compact', 'Padding', 'compact');

    for i = 1:numel(emResults)
        start = emResults(i).emResult.bestStart;

        ax = nexttile(tl, i);
        plot(ax, start.logLikHistory, '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
        title(ax, sprintf('%s | monotone=%d', emResults(i).versionName, ...
            start.logLikelihoodIsMonotone), 'FontWeight', 'normal');
        ylabel(ax, 'log likelihood');
        grid(ax, 'on');

        ax = nexttile(tl, i + numel(emResults));
        semilogy(ax, start.qHistory(:, 1), '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
        hold(ax, 'on');
        semilogy(ax, start.qHistory(:, 2), '-s', 'LineWidth', 1.2, 'MarkerSize', 4);
        ylabel(ax, 'q');
        legend(ax, {'q_{f0}', 'q_{group}'}, 'Location', 'best');
        grid(ax, 'on');

        ax = nexttile(tl, i + 2 * numel(emResults));
        plot(ax, start.rHistory, '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
        xlabel(ax, 'EM iteration');
        ylabel(ax, 'r (dB^2)');
        grid(ax, 'on');
    end

    title(tl, sprintf('%s-band scalar-r EM diagnostics', title_case(bandConfig.bandName)), ...
        'FontWeight', 'normal');
    outFile = fullfile(figuresDir, sprintf('%s_%s_em_diagnostics.png', ...
        bandConfig.outputPrefix, lower(bandConfig.bandName)));
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
end

function plot_profile_q_heatmap(versionResults, bandConfig, figuresDir)
    idx = find([versionResults.versionName] == "Profile-q", 1);
    if isempty(idx) || ~isfield(versionResults(idx), 'profileLogLikSurface') || ...
            isempty(versionResults(idx).profileLogLikSurface)
        return;
    end

    result = versionResults(idx);
    fig = figure('Color', 'w', 'Position', [100 100 780 640]);
    ax = axes(fig);
    imagesc(ax, log10(result.profileQGroupDifferenceGrid), ...
        log10(result.profileQF0Grid), result.profileLogLikSurface);
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    plot(ax, log10(result.qGroupDifference), log10(result.qF0), 'rp', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 16);
    colorbar(ax);
    xlabel(ax, 'log_{10}(q_{group difference})');
    ylabel(ax, 'log_{10}(q_{f0})');
    title(ax, {sprintf('%s-band Profile-q heatmap', title_case(bandConfig.bandName)), ...
        sprintf('r fixed at %.3g dB^2; red star = max', result.profileR)}, ...
        'FontWeight', 'normal');

    outFile = fullfile(figuresDir, sprintf('%s_%s_profile_q_heatmap.png', ...
        bandConfig.outputPrefix, lower(bandConfig.bandName)));
    exportgraphics(fig, outFile, 'Resolution', 200);
    close(fig);
end

function plot_band(ax, ages, meanVals, lowVals, highVals, color)
    ages = ages(:);
    meanVals = meanVals(:);
    lowVals = lowVals(:);
    highVals = highVals(:);
    patch(ax, [ages; flipud(ages)], [lowVals; flipud(highVals)], color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.32);
    lineColor = max(color * 0.55, [0 0 0]);
    plot(ax, ages, meanVals, '-', 'Color', lineColor, 'LineWidth', 2.2);
end

function result = get_version_result(versionResults, versionName)
    idx = find([versionResults.versionName] == string(versionName), 1);
    if isempty(idx)
        error('Version %s was not found.', versionName);
    end
    result = versionResults(idx);
end

function style_axis(ax)
    grid(ax, 'on');
    box(ax, 'off');
    ax.FontName = 'Helvetica';
    ax.FontSize = 11;
    ax.TickDir = 'out';
end

function qInit = compute_q_init(ageYears, powerDB)
    ageSpan = max(ageYears) - min(ageYears);
    powerRange = max(powerDB) - min(powerDB);
    qInit = (powerRange ^ 2) / max(ageSpan ^ 3, eps);
end

function [ageLow, ageHigh] = get_age_range(ageYears, mask)
    if any(mask)
        ageLow = min(ageYears(mask));
        ageHigh = max(ageYears(mask));
    else
        ageLow = NaN;
        ageHigh = NaN;
    end
end

function out = title_case(value)
    value = char(string(value));
    if isempty(value)
        out = value;
    else
        out = [upper(value(1)), lower(value(2:end))];
    end
end
