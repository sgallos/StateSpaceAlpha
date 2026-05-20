%% run_ss_age_diff
% Step 3 entry point for the StateSpaceAlpha analysis.
%
% This script loads the 19-row smoke-test alpha table and fits the 4-D
% baseline + CP-control difference integrated-Wiener-process state-space
% model.
%
% State:
%   x(a) = [baseline; baseline_slope; delta; delta_slope]
%
% Observation:
%   Control: alpha = baseline
%   CP:      alpha = baseline + delta
%
% The primary scientific output is delta(a), the CP - Control alpha-power
% difference over age.

clear; clc;

%% Step 1: Locate the alpha table
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));

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

%% Step 2: Fit the main Step 3 model
fprintf('\nRunning Step 3 baseline/delta smoother with qScaleF0 = 1, qScaleDelta = 1\n');
mainResult = SS_age_diff(alphaTable, ...
    'qScaleF0', 1, ...
    'qScaleDelta', 1, ...
    'verbose', true);

%% Step 3: Run q-sensitivity grid for delta(a)
qScalesF0 = [0.1 1 10];
qScalesDelta = [0.1 1 10];
sensitivityResults = cell(numel(qScalesF0), numel(qScalesDelta));

for i = 1:numel(qScalesF0)
    for j = 1:numel(qScalesDelta)
        fprintf('\nSensitivity fit qScaleF0 = %.3g, qScaleDelta = %.3g\n', ...
            qScalesF0(i), qScalesDelta(j));
        sensitivityResults{i, j} = SS_age_diff(alphaTable, ...
            'qScaleF0', qScalesF0(i), ...
            'qScaleDelta', qScalesDelta(j), ...
            'verbose', false);
    end
end

%% Step 4: Save Step 3 outputs
step3MatFile = fullfile(outDir, 'ssm_step3_age_difference_results.mat');
step3CsvFile = fullfile(outDir, 'ssm_step3_age_difference_trajectory.csv');

writetable(mainResult.trajectory, step3CsvFile);
save(step3MatFile, 'alphaTable', 'mainResult', 'qScalesF0', ...
    'qScalesDelta', 'sensitivityResults', 'alphaTableFileName', 'alphaTableFile');

fprintf('\nSaved Step 3 trajectory CSV to %s\n', step3CsvFile);
fprintf('Saved Step 3 MAT results to %s\n', step3MatFile);

%% Step 5: Figure A, baseline/control and CP trajectories
figOverlay = figure('Color', 'w', 'Position', [100 100 1100 580]);
ax = axes(figOverlay);
hold(ax, 'on');

plot_subject_scatter(ax, mainResult.trajectory);
hControlBand = plot_ci_band(ax, mainResult.trajectory.ageYears, ...
    mainResult.trajectory.baselineCILow_dB, mainResult.trajectory.baselineCIHigh_dB, ...
    [0.72 0.82 1.00], 0.26);
hCPBand = plot_ci_band(ax, mainResult.trajectory.ageYears, ...
    mainResult.trajectory.cpCILow_dB, mainResult.trajectory.cpCIHigh_dB, ...
    [1.00 0.72 0.72], 0.26);
hControlCurve = plot(ax, mainResult.trajectory.ageYears, mainResult.trajectory.baselineMean_dB, ...
    '-', 'Color', [0.05 0.25 0.75], 'LineWidth', 2.5);
hCPCurve = plot(ax, mainResult.trajectory.ageYears, mainResult.trajectory.cpMean_dB, ...
    '-', 'Color', [0.75 0.12 0.12], 'LineWidth', 2.5);

xlabel(ax, 'Age (years)');
ylabel(ax, 'Alpha power (dB)');
title(ax, 'Step 3 baseline/control and CP trajectories');
legend(ax, [hControlBand, hCPBand, hControlCurve, hCPCurve], ...
    {'Control 95% CI', 'CP 95% CI', 'Control / baseline', 'CP = baseline + difference'}, ...
    'Location', 'best');
style_ssm_axes(ax);

overlayPlotFile = fullfile(outDir, 'ssm_step3_baseline_cp_overlay.png');
saveas(figOverlay, overlayPlotFile);
fprintf('Saved Step 3 overlay plot to %s\n', overlayPlotFile);

%% Step 6: Figure B, primary delta(a) scientific plot
figDelta = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figDelta);
hold(ax, 'on');

hDeltaBand = plot_ci_band(ax, mainResult.trajectory.ageYears, ...
    mainResult.trajectory.deltaCILow_dB, mainResult.trajectory.deltaCIHigh_dB, ...
    [0.78 0.80 0.86], 0.40);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hDelta = plot(ax, mainResult.trajectory.ageYears, mainResult.trajectory.deltaMean_dB, ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, 'Primary Step 3 output: CP-Control alpha-power difference across age');
legend(ax, [hDeltaBand, hDelta, hZero], {'Delta 95% CI', 'Delta mean', 'Zero line'}, ...
    'Location', 'best');
style_ssm_axes(ax);

deltaPlotFile = fullfile(outDir, 'ssm_step3_delta_primary.png');
saveas(figDelta, deltaPlotFile);
fprintf('Saved Step 3 primary delta plot to %s\n', deltaPlotFile);

%% Step 7: Figure C, delta(a) q-sensitivity
figSensitivity = figure('Color', 'w', 'Position', [100 100 1120 590]);
ax = axes(figSensitivity);
hold(ax, 'on');

yline(ax, 0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2);

lineHandles = gobjects(numel(qScalesF0) * numel(qScalesDelta), 1);
legendLabels = strings(numel(lineHandles), 1);
lineColors = lines(numel(lineHandles));
lineIndex = 0;

for i = 1:numel(qScalesF0)
    for j = 1:numel(qScalesDelta)
        lineIndex = lineIndex + 1;
        thisResult = sensitivityResults{i, j};
        lineHandles(lineIndex) = plot(ax, thisResult.trajectory.ageYears, ...
            thisResult.trajectory.deltaMean_dB, '-', ...
            'Color', lineColors(lineIndex, :), 'LineWidth', 1.8);
        legendLabels(lineIndex) = sprintf('f0 %.1g, delta %.1g', ...
            qScalesF0(i), qScalesDelta(j));
    end
end

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, 'Step 3 sensitivity: delta(a) across smoothness settings');
legend(ax, lineHandles, cellstr(legendLabels), 'Location', 'eastoutside');
style_ssm_axes(ax);

sensitivityPlotFile = fullfile(outDir, 'ssm_step3_delta_q_sensitivity.png');
saveas(figSensitivity, sensitivityPlotFile);
fprintf('Saved Step 3 delta q-sensitivity plot to %s\n', sensitivityPlotFile);

%% Local plotting helpers
function plot_subject_scatter(ax, trajectory)
    controlMask = strcmpi(trajectory.groupLabel, 'Control');
    cpMask = strcmpi(trajectory.groupLabel, 'CP') | strcmpi(trajectory.groupLabel, 'Patient');

    scatter(ax, trajectory.ageYears(controlMask), trajectory.alpha_dB(controlMask), 70, ...
        'MarkerFaceColor', [0.05 0.25 0.75], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    scatter(ax, trajectory.ageYears(cpMask), trajectory.alpha_dB(cpMask), 70, ...
        'MarkerFaceColor', [0.75 0.12 0.12], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
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
