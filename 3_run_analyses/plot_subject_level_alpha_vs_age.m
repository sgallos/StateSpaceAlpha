%% plot_subject_level_alpha_vs_age
% Plot the subject-level alpha values that enter the state-space model.
%
% Input:
%   outputs/alpha_table_for_ssm_subepochs_collapsed.csv
%
% Output:
%   outputs/subject_level_alpha_vs_age.png

clear; clc;

%% Step 1: Locate inputs and outputs
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));

alphaTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs_collapsed.csv');
plotFile = fullfile(outDir, 'subject_level_alpha_vs_age.png');

if ~isfile(alphaTableFile)
    error(['Missing collapsed alpha table:\n%s\n\n' ...
        'Run 1_data_preparation/collapse_subepoch_alpha_table_for_ssm.m first.'], alphaTableFile);
end

%% Step 2: Load the 19-row subject table
alphaTable = readtable(alphaTableFile, 'TextType', 'string');
alphaTable = sortrows(alphaTable, {'groupLabel','ageYears','subjectID'});

isControl = strcmpi(alphaTable.groupLabel, "Control");
isCP = strcmpi(alphaTable.groupLabel, "CP");

fprintf('Loaded %d subject-level alpha rows from %s\n', height(alphaTable), alphaTableFile);
fprintf('Controls: %d | CP: %d\n', sum(isControl), sum(isCP));
fprintf('Control mean alpha: %.2f dB\n', mean(alphaTable.alpha_dB(isControl), 'omitnan'));
fprintf('CP mean alpha: %.2f dB\n', mean(alphaTable.alpha_dB(isCP), 'omitnan'));

%% Step 3: Plot alpha vs age by group
fig = figure('Color', 'w', 'Position', [100 100 900 560]);
ax = axes(fig);
hold(ax, 'on');

controlColor = [0.12 0.32 0.78];
cpColor = [0.78 0.18 0.16];

scatter(ax, alphaTable.ageYears(isControl), alphaTable.alpha_dB(isControl), ...
    78, 'o', 'MarkerFaceColor', controlColor, 'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.9, 'DisplayName', 'Control');
scatter(ax, alphaTable.ageYears(isCP), alphaTable.alpha_dB(isCP), ...
    78, 'o', 'MarkerFaceColor', cpColor, 'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.9, 'DisplayName', 'CP');

plot_group_mean_line(ax, alphaTable.ageYears(isControl), alphaTable.alpha_dB(isControl), controlColor);
plot_group_mean_line(ax, alphaTable.ageYears(isCP), alphaTable.alpha_dB(isCP), cpColor);

xlabel(ax, 'Age (years)');
ylabel(ax, 'F3 alpha power (dB)');
title(ax, 'Subject-level alpha power vs age at F3');
legend(ax, 'Location', 'best');
style_subject_scatter_axes(ax);

saveas(fig, plotFile);
fprintf('Saved subject-level alpha scatter to %s\n', plotFile);

%% Local helpers
function plot_group_mean_line(ax, ageYears, alphaDB, color)
    if isempty(ageYears)
        return;
    end

    xRange = [min(ageYears), max(ageYears)];
    yMean = mean(alphaDB, 'omitnan');
    plot(ax, xRange, [yMean, yMean], '-', 'Color', color, ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
end

function style_subject_scatter_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
