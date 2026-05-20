%% collapse_subepoch_alpha_table_for_ssm
% Collapse the 76-row 4 x 30 s sub-epoch alpha table into the 19-row table
% format expected by the existing StateSpaceAlpha SSM scripts.
%
% The default observation variance is the standard error of each subject's
% four sub-epoch alpha estimates:
%
%   mfdb_var = var(alpha_subepochs) / nSubEpochs
%
% This intentionally captures within-recording alpha fluctuation, not only
% the much smaller within-sub-epoch MFDB spectral-estimation uncertainty.

clear; clc;

%% Step 1: User-facing parameters
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));

subepochAlphaCsv = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'alpha_table_for_ssm_subepochs.csv');

varianceRule = "empirical_sem";  % "empirical_sem" or "max_empirical_internal"
activeTableCopy = false;         % true overwrites outputs/alpha_table_for_ssm.csv

collapsedCsv = fullfile(outDir, 'alpha_table_for_ssm_subepochs_collapsed.csv');
collapsedMat = fullfile(outDir, 'alpha_table_for_ssm_subepochs_collapsed.mat');
comparisonCsv = fullfile(outDir, 'alpha_table_subepoch_vs_120s_comparison.csv');
comparisonPlot = fullfile(outDir, 'alpha_table_subepoch_vs_120s_comparison.png');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~isfile(subepochAlphaCsv)
    error(['Sub-epoch alpha table not found:\n%s\n\n' ...
        'Run these first in MATLAB_Multitaper_Hz_Domain_BTS:\n' ...
        '  run_subject_mfdb_bootstrap_subepochs_batch\n' ...
        '  smoke_test_alpha_extraction_subepochs'], subepochAlphaCsv);
end

%% Step 2: Collapse 76 sub-epoch rows into 19 subject rows
subepochTable = readtable(subepochAlphaCsv, 'TextType', 'string');
collapsedTable = collapse_subepochs_v2(subepochTable, varianceRule);

writetable(collapsedTable, collapsedCsv);
save(collapsedMat, 'collapsedTable', 'subepochTable', 'varianceRule');

fprintf('Loaded %d sub-epoch rows from %s\n', height(subepochTable), subepochAlphaCsv);
fprintf('Collapsed to %d subject rows using varianceRule = %s\n', ...
    height(collapsedTable), varianceRule);
fprintf('Saved collapsed table to %s\n', collapsedCsv);
disp(collapsedTable(:, {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var', ...
    'withinSubjectAlphaSD_dB','meanSubepochMfdbSD_dB','nSubEpochs'}));

%% Step 3: Optional comparison against the original 120 s table
originalTableFile = fullfile(outDir, 'alpha_table_for_ssm.csv');
comparisonTable = add_original_120s_comparison(collapsedTable, originalTableFile);
writetable(comparisonTable, comparisonCsv);
fprintf('Saved 120 s comparison table to %s\n', comparisonCsv);

if isfile(originalTableFile)
    make_comparison_plot(comparisonTable, comparisonPlot);
    fprintf('Saved comparison plot to %s\n', comparisonPlot);
else
    fprintf('Original 120 s alpha table not found; skipped comparison plot.\n');
end

%% Step 4: Optional active-table copy
if activeTableCopy
    activeCsv = fullfile(outDir, 'alpha_table_for_ssm.csv');
    writetable(collapsedTable, activeCsv);
    fprintf('Copied collapsed table to active SSM table: %s\n', activeCsv);
else
    fprintf(['Active SSM table was not overwritten.\n' ...
        'To run Step 3/4 on this table, set alphaTableFileName to:\n' ...
        '  alpha_table_for_ssm_subepochs_collapsed.csv\n']);
end

%% Local functions
function T = collapse_subepochs_v2(subepochTable, varianceRule)
    requiredVars = {'subjectID','subepochIdx','groupLabel','ageYears','alpha_dB','mfdb_var'};
    missingVars = setdiff(requiredVars, subepochTable.Properties.VariableNames);
    if ~isempty(missingVars)
        error('Sub-epoch table is missing required columns: %s', strjoin(missingVars, ', '));
    end

    subepochTable = sortrows(subepochTable, {'groupLabel','ageYears','subjectID','subepochIdx'});
    subjectIDs = unique(subepochTable.subjectID, 'stable');
    nSubjects = numel(subjectIDs);

    T = table('Size', [nSubjects, 15], ...
        'VariableTypes', {'string','string','double','double','double','double','double', ...
        'double','double','double','double','double','double','double','string'}, ...
        'VariableNames', {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var', ...
        'nSubEpochs','alphaSimpleMean_dB','alphaWeightedMean_dB', ...
        'withinSubjectAlphaSD_dB','withinSubjectAlphaRange_dB', ...
        'empiricalSemVariance','meanSubepochMfdbVariance', ...
        'weightedInternalSemVariance','meanSubepochMfdbSD_dB','varianceRule'});

    for i = 1:nSubjects
        rows = subepochTable.subjectID == subjectIDs(i);
        thisData = sortrows(subepochTable(rows, :), 'subepochIdx');

        alphaVals = thisData.alpha_dB(:);
        subepochVars = max(thisData.mfdb_var(:), eps);
        nSubEpochs = numel(alphaVals);
        weights = 1 ./ subepochVars;

        alphaSimpleMean = mean(alphaVals, 'omitnan');
        alphaWeightedMean = sum(weights .* alphaVals) / sum(weights);
        withinSubjectSD = std(alphaVals, 0, 'omitnan');
        withinSubjectRange = max(alphaVals) - min(alphaVals);
        empiricalSemVariance = var(alphaVals, 0, 'omitnan') / nSubEpochs;
        meanSubepochMfdbVariance = mean(subepochVars, 'omitnan');
        weightedInternalSemVariance = 1 / sum(weights);

        switch lower(string(varianceRule))
            case "empirical_sem"
                effectiveVariance = empiricalSemVariance;

            case "max_empirical_internal"
                effectiveVariance = max(empiricalSemVariance, weightedInternalSemVariance);

            otherwise
                error('Unknown varianceRule "%s".', varianceRule);
        end

        T.subjectID(i) = subjectIDs(i);
        T.groupLabel(i) = thisData.groupLabel(1);
        T.ageYears(i) = thisData.ageYears(1);
        T.alpha_dB(i) = alphaSimpleMean;
        T.mfdb_var(i) = max(effectiveVariance, eps);
        T.nSubEpochs(i) = nSubEpochs;
        T.alphaSimpleMean_dB(i) = alphaSimpleMean;
        T.alphaWeightedMean_dB(i) = alphaWeightedMean;
        T.withinSubjectAlphaSD_dB(i) = withinSubjectSD;
        T.withinSubjectAlphaRange_dB(i) = withinSubjectRange;
        T.empiricalSemVariance(i) = empiricalSemVariance;
        T.meanSubepochMfdbVariance(i) = meanSubepochMfdbVariance;
        T.weightedInternalSemVariance(i) = weightedInternalSemVariance;
        T.meanSubepochMfdbSD_dB(i) = mean(sqrt(subepochVars), 'omitnan');
        T.varianceRule(i) = string(varianceRule);
    end
end

function comparisonTable = add_original_120s_comparison(collapsedTable, originalTableFile)
    comparisonTable = collapsedTable;
    comparisonTable.original120sAlpha_dB = nan(height(comparisonTable), 1);
    comparisonTable.original120sMfdbVar = nan(height(comparisonTable), 1);
    comparisonTable.alphaDiffCollapsedMinus120s_dB = nan(height(comparisonTable), 1);
    comparisonTable.mfdbVarRatioCollapsedOver120s = nan(height(comparisonTable), 1);

    if ~isfile(originalTableFile)
        return;
    end

    originalTable = readtable(originalTableFile, 'TextType', 'string');
    for i = 1:height(comparisonTable)
        matchIdx = find(strcmpi(originalTable.subjectID, comparisonTable.subjectID(i)), 1, 'first');
        if isempty(matchIdx)
            continue;
        end

        comparisonTable.original120sAlpha_dB(i) = originalTable.alpha_dB(matchIdx);
        comparisonTable.original120sMfdbVar(i) = originalTable.mfdb_var(matchIdx);
        comparisonTable.alphaDiffCollapsedMinus120s_dB(i) = ...
            comparisonTable.alpha_dB(i) - originalTable.alpha_dB(matchIdx);
        comparisonTable.mfdbVarRatioCollapsedOver120s(i) = ...
            comparisonTable.mfdb_var(i) / originalTable.mfdb_var(matchIdx);
    end
end

function make_comparison_plot(comparisonTable, plotFile)
    hasOriginal = isfinite(comparisonTable.original120sAlpha_dB);
    if ~any(hasOriginal)
        return;
    end

    fig = figure('Color', 'w', 'Position', [100 100 1150 500]);
    tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    controlMask = strcmpi(comparisonTable.groupLabel, 'Control');
    cpMask = strcmpi(comparisonTable.groupLabel, 'CP') | strcmpi(comparisonTable.groupLabel, 'Patient');
    controlColor = [0.9 0.45 0.1];
    cpColor = [0.2 0.45 0.85];

    ax1 = nexttile;
    hold(ax1, 'on');
    plot_group_pair(ax1, comparisonTable, controlMask, controlColor);
    plot_group_pair(ax1, comparisonTable, cpMask, cpColor);
    xlabel(ax1, 'Age (years)');
    ylabel(ax1, 'Alpha power (dB)');
    title(ax1, '120 s vs collapsed 4 x 30 s alpha');
    legend(ax1, {'Control', 'CP'}, 'Location', 'best');
    style_axis(ax1);

    ax2 = nexttile;
    hold(ax2, 'on');
    scatter(ax2, comparisonTable.ageYears(controlMask), ...
        comparisonTable.mfdbVarRatioCollapsedOver120s(controlMask), 70, ...
        'MarkerFaceColor', controlColor, 'MarkerEdgeColor', 'k');
    scatter(ax2, comparisonTable.ageYears(cpMask), ...
        comparisonTable.mfdbVarRatioCollapsedOver120s(cpMask), 70, ...
        'MarkerFaceColor', cpColor, 'MarkerEdgeColor', 'k');
    yline(ax2, 1, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2);
    xlabel(ax2, 'Age (years)');
    ylabel(ax2, 'Collapsed mfdb\_var / 120 s mfdb\_var');
    title(ax2, 'Effective variance change');
    style_axis(ax2);

    saveas(fig, plotFile);
end

function plot_group_pair(ax, T, mask, color)
    idx = find(mask & isfinite(T.original120sAlpha_dB));
    for j = 1:numel(idx)
        i = idx(j);
        plot(ax, [T.ageYears(i) T.ageYears(i)], ...
            [T.original120sAlpha_dB(i) T.alpha_dB(i)], '-', ...
            'Color', [color 0.35], 'HandleVisibility', 'off');
    end
    scatter(ax, T.ageYears(idx), T.original120sAlpha_dB(idx), 90, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', color, ...
        'LineWidth', 1.6, 'HandleVisibility', 'off');
    scatter(ax, T.ageYears(idx), T.alpha_dB(idx), 46, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5);
end

function style_axis(ax)
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'FontName', 'Helvetica', 'FontSize', 11, ...
        'LineWidth', 1, 'TickDir', 'out');
end
