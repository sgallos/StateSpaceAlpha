%% check_subepoch_consistency
% Summarize whether each subject's four 30 s alpha estimates are stable.
%
% Input:
%   outputs/alpha_table_for_ssm_subepochs.csv
%
% Outputs:
%   outputs/subepoch_consistency_summary.csv
%   outputs/subepoch_consistency_scatter.png

clear; clc;

%% Step 1: Locate inputs and outputs
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');

localSubepochTable = fullfile(outDir, 'alpha_table_for_ssm_subepochs.csv');
upstreamSubepochTable = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'alpha_table_for_ssm_subepochs.csv');

summaryFile = fullfile(outDir, 'subepoch_consistency_summary.csv');
plotFile = fullfile(outDir, 'subepoch_consistency_scatter.png');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~isfile(localSubepochTable)
    if isfile(upstreamSubepochTable)
        copyfile(upstreamSubepochTable, localSubepochTable);
        fprintf('Copied sub-epoch table into StateSpaceAlpha outputs:\n  %s\n', localSubepochTable);
    else
        error(['Sub-epoch alpha table not found:\n%s\n\n' ...
            'Run smoke_test_alpha_extraction_subepochs in MATLAB_Multitaper_Hz_Domain_BTS first.'], ...
            localSubepochTable);
    end
end

%% Step 2: Compute per-subject consistency summaries
T = readtable(localSubepochTable, 'TextType', 'string');
requiredVars = {'subjectID','subepochIdx','groupLabel','ageYears','alpha_dB','mfdb_var'};
missingVars = setdiff(requiredVars, T.Properties.VariableNames);
if ~isempty(missingVars)
    error('Sub-epoch table is missing required columns: %s', strjoin(missingVars, ', '));
end

subjectIDs = unique(T.subjectID, 'stable');
nSubjects = numel(subjectIDs);

perSubject = table('Size', [nSubjects, 8], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double'}, ...
    'VariableNames', {'subjectID','groupLabel','ageYears','nSubEpochs', ...
    'alphaMean_dB','alphaSD_dB','alphaRange_dB','meanMfdbSD_dB'});

for i = 1:nSubjects
    rows = T.subjectID == subjectIDs(i);
    thisData = sortrows(T(rows, :), 'subepochIdx');
    alphaVals = thisData.alpha_dB(:);

    perSubject.subjectID(i) = subjectIDs(i);
    perSubject.groupLabel(i) = thisData.groupLabel(1);
    perSubject.ageYears(i) = thisData.ageYears(1);
    perSubject.nSubEpochs(i) = height(thisData);
    perSubject.alphaMean_dB(i) = mean(alphaVals, 'omitnan');
    perSubject.alphaSD_dB(i) = std(alphaVals, 0, 'omitnan');
    perSubject.alphaRange_dB(i) = max(alphaVals) - min(alphaVals);
    perSubject.meanMfdbSD_dB(i) = mean(sqrt(thisData.mfdb_var), 'omitnan');
end

perSubject = sortrows(perSubject, 'alphaSD_dB', 'descend');
writetable(perSubject, summaryFile);

fprintf('Per-subject sub-epoch consistency:\n');
disp(perSubject);
fprintf('\nMedian within-subject SD: %.3f dB\n', median(perSubject.alphaSD_dB, 'omitnan'));
fprintf('Median within-subject range: %.3f dB\n', median(perSubject.alphaRange_dB, 'omitnan'));
fprintf('Median mean MFDB SD: %.4f dB\n', median(perSubject.meanMfdbSD_dB, 'omitnan'));
fprintf('Saved consistency summary to %s\n', summaryFile);

%% Step 3: Plot the four sub-epoch values for each subject
fig = figure('Color', 'w', 'Position', [100 100 1050 620]);
ax = axes(fig);
hold(ax, 'on');

controlColor = [0.12 0.32 0.78];
cpColor = [0.78 0.18 0.16];

for i = 1:nSubjects
    rows = T.subjectID == subjectIDs(i);
    thisData = sortrows(T(rows, :), 'subepochIdx');

    if strcmpi(thisData.groupLabel(1), "CP")
        thisColor = cpColor;
    else
        thisColor = controlColor;
    end

    jitter = (thisData.subepochIdx - mean(thisData.subepochIdx)) * 0.045;
    plot(ax, thisData.ageYears + jitter, thisData.alpha_dB, 'o-', ...
        'Color', thisColor, 'LineWidth', 0.9, 'MarkerSize', 4, ...
        'HandleVisibility', 'off');
end

hControl = plot(ax, nan, nan, 'o-', 'Color', controlColor, ...
    'LineWidth', 1.2, 'MarkerSize', 5);
hCP = plot(ax, nan, nan, 'o-', 'Color', cpColor, ...
    'LineWidth', 1.2, 'MarkerSize', 5);

xlabel(ax, 'Age (years), jittered by sub-epoch');
ylabel(ax, 'Alpha power (dB)');
title(ax, 'Sub-epoch alpha estimates per subject');
legend(ax, [hControl hCP], {'Control', 'CP'}, 'Location', 'best');
style_consistency_axes(ax);

saveas(fig, plotFile);
fprintf('Saved consistency scatter plot to %s\n', plotFile);

%% Local helpers
function style_consistency_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
