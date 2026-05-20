function plot_subjects_with_trajectories(resultsFile, outputFile)
%PLOT_SUBJECTS_WITH_TRAJECTORIES Plot data with Control and CP SSM curves.
%
% Shows the 19 subject-level alpha values overlaid on the posterior Control
% trajectory f0(a) and the posterior CP trajectory f0(a) + delta(a).

if nargin < 1 || isempty(resultsFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    resultsFile = fullfile(repoRoot, 'outputs', 'ssm_posterior_no_resampling_results.mat');
end

if nargin < 2 || isempty(outputFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    outputFile = fullfile(repoRoot, 'outputs', 'subjects_with_trajectories.png');
end

if ~isfile(resultsFile)
    error(['Missing primary SSM results file:\n%s\n\n' ...
        'Run 3_run_analyses/run_ssm_posterior_no_resampling.m first.'], resultsFile);
end

S = load(resultsFile);
if ~isfield(S, 'primaryResult') || ~isfield(S.primaryResult, 'trajectory')
    error('Expected results file to contain primaryResult.trajectory.');
end

trajectory = sortrows(S.primaryResult.trajectory, {'ageYears','subjectID'});
assert_required_columns(trajectory);

ageYears = trajectory.ageYears(:);
alphaDB = trajectory.alpha_dB(:);
isCP = strcmpi(trajectory.groupLabel, "CP");
isControl = ~isCP;

controlMean = trajectory.baselineMean_dB(:);
controlLow = trajectory.baselineCILow_dB(:);
controlHigh = trajectory.baselineCIHigh_dB(:);

cpMean = trajectory.cpMean_dB(:);
cpLow = trajectory.cpCILow_dB(:);
cpHigh = trajectory.cpCIHigh_dB(:);

controlColor = [0.12 0.32 0.78];
cpColor = [0.78 0.18 0.16];
controlBandColor = [0.68 0.76 0.92];
cpBandColor = [0.93 0.72 0.70];

fig = figure('Color', 'w', 'Position', [100 100 1050 620]);
ax = axes(fig);
hold(ax, 'on');

hControlBand = plot_ci_band(ax, ageYears, controlLow, controlHigh, controlBandColor, 0.38);
hCPBand = plot_ci_band(ax, ageYears, cpLow, cpHigh, cpBandColor, 0.38);

hControlLine = plot(ax, ageYears, controlMean, '-', 'Color', controlColor, ...
    'LineWidth', 2.4);
hCPLine = plot(ax, ageYears, cpMean, '-', 'Color', cpColor, ...
    'LineWidth', 2.4);

hControlSubjects = scatter(ax, ageYears(isControl), alphaDB(isControl), ...
    65, controlColor, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);
hCPSubjects = scatter(ax, ageYears(isCP), alphaDB(isCP), ...
    65, cpColor, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);

xlabel(ax, 'Age (years)');
ylabel(ax, 'F3 alpha power (dB)');
title(ax, {'Control and CP alpha trajectories with subject-level data', ...
    'Bands are 95% posterior credible intervals from the SSM smoother'}, ...
    'FontWeight', 'normal');
legend(ax, [hControlBand, hCPBand, hControlLine, hCPLine, hControlSubjects, hCPSubjects], ...
    {'Control 95% credible band', 'CP 95% credible band', ...
    'Control trajectory f0(a)', 'CP trajectory f0(a) + delta(a)', ...
    'Control subjects', 'CP subjects'}, ...
    'Location', 'best', 'Interpreter', 'none');
style_subject_trajectory_axes(ax);

outputDir = fileparts(outputFile);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
exportgraphics(fig, outputFile, 'Resolution', 200);
fprintf('Saved subjects-with-trajectories figure to %s\n', outputFile);
end

function assert_required_columns(trajectory)
    requiredColumns = {'subjectID','groupLabel','ageYears','alpha_dB', ...
        'observationVariance','baselineMean_dB','baselineCILow_dB', ...
        'baselineCIHigh_dB','cpMean_dB','cpCILow_dB','cpCIHigh_dB'};
    missingColumns = setdiff(requiredColumns, trajectory.Properties.VariableNames);
    if ~isempty(missingColumns)
        error('Trajectory table is missing required columns: %s', strjoin(missingColumns, ', '));
    end
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

function style_subject_trajectory_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
