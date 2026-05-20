%% run_ss_age_diff
% Step 2 entry point for the StateSpaceAlpha analysis.
%
% This script loads the 19-row smoke-test alpha table and fits a pooled
% 2-D integrated-Wiener-process state-space smoother. CP/control labels are
% plotted for visual context only; they are not used by the model yet.

clear; clc;

%% Step 1: Locate the smoke-test alpha table
repoRoot = fileparts(mfilename('fullpath'));
outDir = fullfile(repoRoot, 'outputs');
alphaTableFile = fullfile(outDir, 'alpha_table_for_ssm.csv');

if ~isfile(alphaTableFile)
    error(['Alpha table not found:\n%s\n\n' ...
        'Run this first:\n' ...
        '  cd(''%s'')\n' ...
        '  smoke_test_alpha_extraction'], alphaTableFile, repoRoot);
end

alphaTable = readtable(alphaTableFile, 'TextType', 'string');
fprintf('Loaded alpha table: %s\n', alphaTableFile);
disp(alphaTable(:, {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var'}));

%% Step 2: Fit q-sensitivity versions of the pooled 2-D IWP model
qScales = [0.1 1 10];
ssmResults = cell(numel(qScales), 1);

for i = 1:numel(qScales)
    fprintf('\nRunning pooled Step 2 smoother with qScale = %.3g\n', qScales(i));
    ssmResults{i} = SS_age_diff(alphaTable, 'qScale', qScales(i), 'verbose', true);
end

mainResult = ssmResults{2}; % qScale = 1

%% Step 3: Save outputs for later steps
step2MatFile = fullfile(outDir, 'ssm_step2_pooled_iwp_results.mat');
step2CsvFile = fullfile(outDir, 'ssm_step2_pooled_iwp_q1_trajectory.csv');

writetable(mainResult.trajectory, step2CsvFile);
save(step2MatFile, 'alphaTable', 'qScales', 'ssmResults', 'mainResult');

fprintf('\nSaved Step 2 trajectory CSV to %s\n', step2CsvFile);
fprintf('Saved Step 2 MAT results to %s\n', step2MatFile);

%% Step 4: Plot the main qScale = 1 fit
figMain = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figMain);
hold(ax, 'on');

plot_subject_scatter(ax, mainResult.trajectory);
hBand = plot_smoother_band(ax, mainResult.trajectory.ageYears, ...
    mainResult.trajectory.smoothCILow_dB, mainResult.trajectory.smoothCIHigh_dB, ...
    [0.72 0.78 0.86]);
hFilter = plot(ax, mainResult.trajectory.ageYears, mainResult.trajectory.filterMean_dB, ...
    '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5);
hSmooth = plot(ax, mainResult.trajectory.ageYears, mainResult.trajectory.smoothMean_dB, ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.5);

xlabel(ax, 'Age (years)');
ylabel(ax, 'Alpha power (dB)');
title(ax, sprintf('Step 2 pooled 2-D IWP smoother, qScale = 1, q = %.3g', ...
    mainResult.diagnostics.processNoiseQ));
legend(ax, [hBand, hFilter, hSmooth], ...
    {'RTS smoother 95% CI', 'Forward filter mean', 'RTS smoother mean'}, ...
    'Location', 'best');
style_ssm_axes(ax);

mainPlotFile = fullfile(outDir, 'ssm_step2_pooled_iwp_q1.png');
saveas(figMain, mainPlotFile);
fprintf('Saved Step 2 main plot to %s\n', mainPlotFile);

%% Step 5: Plot q-sensitivity
figSensitivity = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figSensitivity);
hold(ax, 'on');

plot_subject_scatter(ax, mainResult.trajectory);

lineStyles = {'--', '-', ':'};
lineColors = [0.45 0.45 0.45; 0.05 0.05 0.05; 0.55 0.15 0.15];
lineHandles = gobjects(numel(qScales), 1);
legendLabels = strings(numel(qScales), 1);

for i = 1:numel(qScales)
    thisTrajectory = ssmResults{i}.trajectory;
    lineHandles(i) = plot(ax, thisTrajectory.ageYears, thisTrajectory.smoothMean_dB, ...
        lineStyles{i}, 'Color', lineColors(i, :), 'LineWidth', 2.3);
    legendLabels(i) = sprintf('qScale %.1g', qScales(i));
end

xlabel(ax, 'Age (years)');
ylabel(ax, 'Alpha power (dB)');
title(ax, 'Step 2 q-sensitivity: pooled 2-D IWP smoother');
legend(ax, lineHandles, cellstr(legendLabels), 'Location', 'best');
style_ssm_axes(ax);

sensitivityPlotFile = fullfile(outDir, 'ssm_step2_q_sensitivity.png');
saveas(figSensitivity, sensitivityPlotFile);
fprintf('Saved Step 2 q-sensitivity plot to %s\n', sensitivityPlotFile);

%% Local plotting helpers
function plot_subject_scatter(ax, trajectory)
    controlMask = strcmpi(trajectory.groupLabel, 'Control');
    cpMask = strcmpi(trajectory.groupLabel, 'CP') | strcmpi(trajectory.groupLabel, 'Patient');

    scatter(ax, trajectory.ageYears(controlMask), trajectory.alpha_dB(controlMask), 70, ...
        'MarkerFaceColor', [0.9 0.45 0.1], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    scatter(ax, trajectory.ageYears(cpMask), trajectory.alpha_dB(cpMask), 70, ...
        'MarkerFaceColor', [0.2 0.45 0.85], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
end

function hBand = plot_smoother_band(ax, xVals, ciLow, ciHigh, bandColor)
    xVals = xVals(:);
    ciLow = ciLow(:);
    ciHigh = ciHigh(:);
    validMask = isfinite(xVals) & isfinite(ciLow) & isfinite(ciHigh);

    hBand = patch(ax, ...
        [xVals(validMask); flipud(xVals(validMask))], ...
        [ciLow(validMask); flipud(ciHigh(validMask))], ...
        bandColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35);
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
