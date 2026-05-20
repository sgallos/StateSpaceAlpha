%% plot_mfdb_alpha_power_vs_age
% Plot MFDB-derived alpha-band power vs age using validated subject-level MFDB outputs.
% 6A = raw subject-level points.
% 6B = age-window alpha power with MFDB-based 95% confidence intervals.
%
% This script supports MFDB-native age-axis analysis. Stationary
% non-overlapping bins are the cleaner inferential mode because each age
% point uses a distinct age interval. Sliding overlapping windows are the
% smoother trajectory-visualization mode. Both modes propagate upstream
% subject-level MFDB uncertainty into the within-window group summary.

%% Step 1: Setup
clear; clc;

cfg = mfdb_config();
addpath(cfg.scriptRoot);
addpath(cfg.pedsRoot);
addpath(fullfile(cfg.pedsRoot, 'src'));
addpath(genpath(fullfile(cfg.pedsRoot, 'src')));

anesthesia = cfg.anesthesia;
electrodeLabel = cfg.electrodeLabel;
bandName = 'Alpha';
bandRange = [8 13];
ageMin = 10;
ageMax = 15.5;

%% Figure settings that change the sensitivity-analysis output
% Change these values to regenerate Figure 3 under different assumptions.
groupSummaryStatistic = 'mean';     % 'median' or 'mean'
stationary = 1.25;
ageWindowHalfWidthYears = 2;       % stationary: 0.5 or 1; sliding: 1.25 or 2
ageWindowStepYears = .5;              % ignored when useStationaryBins is true

% Windowing mode for the age axis. Both modes propagate MFDB uncertainty
% into the within-window group summary.
%   true  = non-overlapping bins: independent pointwise CIs, smaller per-bin n
%   false = sliding overlapping windows: smoothed trajectory, pointwise CIs only
useStationaryBins = false;

validGroupSummaryStatistics = {'median', 'mean'};
if ~any(strcmpi(groupSummaryStatistic, validGroupSummaryStatistics))
    error('groupSummaryStatistic must be ''median'' or ''mean''.');
end
groupSummaryStatistic = lower(groupSummaryStatistic);
groupSummaryFunction = str2func(groupSummaryStatistic);

if useStationaryBins
    ageWindowModeTag = 'stationary';
    ageWindowModeLabel = 'Stationary Bin';
else
    ageWindowModeTag = 'sliding';
    ageWindowModeLabel = 'Sliding Window';
end

outDir = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_alpha_power_vs_age');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

plotFileBase = sprintf('mfdb_%s_%s_%0.1f_to_%0.1f_%s_%s_%s_hw%0.2f_step%0.2f', ...
    lower(bandName), lower(electrodeLabel), ageMin, ageMax, anesthesia, ...
    lower(groupSummaryStatistic), ageWindowModeTag, ageWindowHalfWidthYears, ageWindowStepYears);
validationFile = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_validation', 'mfdb_validation_results.mat');

fprintf('Using validated MFDB subject-level outputs\n');
fprintf('Electrode request: %s\n', electrodeLabel);
fprintf('Band: %s [%g %g] Hz\n', bandName, bandRange(1), bandRange(2));
fprintf('Age range: [%.1f %.1f]\n', ageMin, ageMax);
fprintf('Group summary statistic: %s\n', groupSummaryStatistic);
fprintf('Age window mode: %s\n', ageWindowModeTag);
fprintf('Age window half-width: %.2f years\n', ageWindowHalfWidthYears);
fprintf('Age window step: %.2f years\n', ageWindowStepYears);
fprintf('Validation file: %s\n', validationFile);

%% Step 2: Load validated subject registry and filter by age
if ~isfile(validationFile)
    error('Validation file not found: %s. Run validate_group_mfdb_inputs.m first.', validationFile);
end

V = load(validationFile, 'includedAfterValidation', 'registry');
if ~isfield(V, 'includedAfterValidation') || isempty(V.includedAfterValidation)
    error('No included subjects found in %s', validationFile);
end

registry = V.includedAfterValidation;
ageMask = isfinite(registry.ageYears) & registry.ageYears >= ageMin & registry.ageYears <= ageMax;
registry = sortrows(registry(ageMask, :), {'groupLabel', 'ageYears', 'subjectID'});

if isempty(registry)
    error('No validated subjects found in age range [%.1f %.1f].', ageMin, ageMax);
end

fprintf('Validated subjects in age range: %d\n', height(registry));
disp(registry(:, {'subjectID', 'groupLabel', 'ageYears', 'usedElectrodeLabel'}));

%% Step 3: Load MFDB subject files and compute subject-level alpha power
results = table();
results.SubjectID = strings(0, 1);
results.Group = strings(0, 1);
results.AgeYears = zeros(0, 1);
results.RequestedElectrode = strings(0, 1);
results.UsedElectrode = strings(0, 1);
results.UsedFallbackElectrode = false(0, 1);
results.BandName = strings(0, 1);
results.BandLowHz = zeros(0, 1);
results.BandHighHz = zeros(0, 1);
results.AlphaPowerLinear = zeros(0, 1);
results.AlphaPower_dB = zeros(0, 1);
results.NBoot = zeros(0, 1);
results.SourceFile = strings(0, 1);

subjectBootPowerDB = cell(height(registry), 1);
subjectBootPowerLinear = cell(height(registry), 1);
commonFreq = [];

for i = 1:height(registry)
    matPath = char(registry.filePath(i));
    S = load(matPath, 'freq', 'S_subject_original_linear', 'S_subject_boot_linear', ...
        'subjectID', 'groupLabel', 'ageYears', 'requestedElectrodeLabel', ...
        'usedElectrodeLabel', 'usedFallbackElectrode');

    freq = S.freq(:);
    if isempty(commonFreq)
        commonFreq = freq;
    elseif ~isequal(size(freq), size(commonFreq)) || any(abs(freq(:) - commonFreq(:)) > 1e-10)
        error('Frequency vector mismatch in %s. Re-run validation.', matPath);
    end

    bandMask = freq >= bandRange(1) & freq <= bandRange(2);
    if ~any(bandMask)
        error('No %s bins found in [%g %g] Hz for %s', bandName, bandRange(1), bandRange(2), matPath);
    end

    originalLinear = trapz(freq(bandMask), S.S_subject_original_linear(bandMask));
    originalLinear = max(originalLinear, eps);
    originalDB = 10 * log10(originalLinear);

    bootLinear = trapz(freq(bandMask), S.S_subject_boot_linear(:, bandMask), 2);
    bootLinear = max(bootLinear, eps);
    bootDB = 10 * log10(bootLinear);

    subjectBootPowerLinear{i} = bootLinear(:);
    subjectBootPowerDB{i} = bootDB(:);

    row = table( ...
        string(S.subjectID), ...
        standardize_group_label(string(S.groupLabel)), ...
        double(S.ageYears), ...
        string(S.requestedElectrodeLabel), ...
        string(S.usedElectrodeLabel), ...
        logical(S.usedFallbackElectrode), ...
        string(bandName), ...
        bandRange(1), ...
        bandRange(2), ...
        originalLinear, ...
        originalDB, ...
        numel(bootDB), ...
        string(matPath), ...
        'VariableNames', {'SubjectID','Group','AgeYears','RequestedElectrode', ...
        'UsedElectrode','UsedFallbackElectrode','BandName','BandLowHz','BandHighHz', ...
        'AlphaPowerLinear','AlphaPower_dB','NBoot','SourceFile'});

    results = [results; row]; %#ok<AGROW>
end

results = sortrows(results, {'Group', 'AgeYears', 'SubjectID'});

fprintf('\nSubject-level MFDB alpha power results:\n');
disp(results(:, {'SubjectID', 'Group', 'AgeYears', 'UsedElectrode', 'AlphaPower_dB'}));

%% Step 4: Build age windows and group-specific subject masks
if useStationaryBins
    ageWindowWidthYears = 2 * ageWindowHalfWidthYears;
    ageWindowCenters = (ageMin + ageWindowHalfWidthYears : ...
        ageWindowWidthYears : ageMax).';
    ageWindowLow = ageWindowCenters - ageWindowHalfWidthYears;
    ageWindowHigh = ageWindowCenters + ageWindowHalfWidthYears;
    fprintf('Mode: STATIONARY non-overlapping bins (bin width = %.2f yr; ageWindowStepYears ignored).\n', ...
        ageWindowWidthYears);
else
    ageWindowWidthYears = 2 * ageWindowHalfWidthYears;
    ageWindowCenters = (ageMin:ageWindowStepYears:ageMax).';
    ageWindowLow = ageWindowCenters - ageWindowHalfWidthYears;
    ageWindowHigh = ageWindowCenters + ageWindowHalfWidthYears;
    fprintf('Mode: SLIDING overlapping windows (half-width = %.2f yr, step = %.2f yr).\n', ...
        ageWindowHalfWidthYears, ageWindowStepYears);
end

controlMask = strcmpi(results.Group, 'Control');
patientMask = strcmpi(results.Group, 'Patient');

controlResults = results(controlMask, :);
patientResults = results(patientMask, :);
controlBoot = subjectBootPowerDB(controlMask);
patientBoot = subjectBootPowerDB(patientMask);

ctrlColor = [0.9 0.45 0.1];
patColor = [0.2 0.45 0.85];
ctrlCIColor = [1.0 0.77 0.58];
patCIColor = [0.64 0.80 1.0];

%% Step 5: Compute sliding-window summaries and MFDB-based confidence intervals
controlWindowSummary = compute_age_window_summary( ...
    controlResults.AgeYears, controlResults.AlphaPower_dB, ...
    ageWindowCenters, ageWindowHalfWidthYears, groupSummaryFunction);
patientWindowSummary = compute_age_window_summary( ...
    patientResults.AgeYears, patientResults.AlphaPower_dB, ...
    ageWindowCenters, ageWindowHalfWidthYears, groupSummaryFunction);

[controlWindowCILow, controlWindowCIHigh, controlWindowSubjectCount] = compute_age_window_mfdb_ci( ...
    controlResults.AgeYears, controlBoot, ageWindowCenters, ...
    ageWindowHalfWidthYears, cfg.Bgroup, groupSummaryFunction);
[patientWindowCILow, patientWindowCIHigh, patientWindowSubjectCount] = compute_age_window_mfdb_ci( ...
    patientResults.AgeYears, patientBoot, ageWindowCenters, ...
    ageWindowHalfWidthYears, cfg.Bgroup, groupSummaryFunction);

windowSupportTable = table( ...
    ageWindowCenters(:), ...
    ageWindowLow(:), ...
    ageWindowHigh(:), ...
    repmat(string(groupSummaryStatistic), numel(ageWindowCenters), 1), ...
    repmat(string(ageWindowModeTag), numel(ageWindowCenters), 1), ...
    repmat(useStationaryBins, numel(ageWindowCenters), 1), ...
    repmat(ageWindowHalfWidthYears, numel(ageWindowCenters), 1), ...
    repmat(ageWindowWidthYears, numel(ageWindowCenters), 1), ...
    repmat(ageWindowStepYears, numel(ageWindowCenters), 1), ...
    controlWindowSubjectCount(:), ...
    patientWindowSubjectCount(:), ...
    controlWindowSummary(:), ...
    patientWindowSummary(:), ...
    controlWindowCILow(:), ...
    controlWindowCIHigh(:), ...
    patientWindowCILow(:), ...
    patientWindowCIHigh(:), ...
    'VariableNames', {'AgeWindowCenter','AgeWindowLow','AgeWindowHigh', ...
    'GroupSummaryStatistic','AgeWindowMode','UseStationaryBins', ...
    'AgeWindowHalfWidthYears','AgeWindowWidthYears','AgeWindowStepYears', ...
    'ControlSubjectCount','PatientSubjectCount', ...
    'ControlSummary_dB','PatientSummary_dB', ...
    'ControlCILow_dB','ControlCIHigh_dB', ...
    'PatientCILow_dB','PatientCIHigh_dB'});

%% Step 6A: Plot raw subject-level alpha power points
figRaw = figure('Color', 'w', 'Position', [100 100 1150 520]);
ax = axes(figRaw);
hold(ax, 'on');

for x = ageMin:1:ceil(ageMax)
    xline(ax, x, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
end

hCtrl = scatter(ax, controlResults.AgeYears, controlResults.AlphaPower_dB, 70, ...
    'MarkerFaceColor', ctrlColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
hPat = scatter(ax, patientResults.AgeYears, patientResults.AlphaPower_dB, 70, ...
    'MarkerFaceColor', patColor, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);

style_axes(ax, ageMin, ageMax);
xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('%s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, sprintf('6A: Raw MFDB %s Power vs Age at %s (%s)', ...
    bandName, upper(electrodeLabel), upper(anesthesia)), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hCtrl, hPat], {'Controls', 'Patients'}, 'Location', 'best');

%% Step 6B: Plot age-window alpha power with MFDB-based 95% CI
figSlidingWindow = figure('Color', 'w', 'Position', [100 100 1150 520]);
ax = axes(figSlidingWindow);
hold(ax, 'on');

for x = ageMin:1:ceil(ageMax)
    xline(ax, x, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
end

hCtrlCI = plot_valid_window_ci(ax, ageWindowCenters, controlWindowSummary, ...
    controlWindowCILow, controlWindowCIHigh, ctrlCIColor);
hPatCI = plot_valid_window_ci(ax, ageWindowCenters, patientWindowSummary, ...
    patientWindowCILow, patientWindowCIHigh, patCIColor);
hCtrlWindow = plot_valid_windows(ax, ageWindowCenters, controlWindowSummary, '-', ctrlColor);
hPatWindow = plot_valid_windows(ax, ageWindowCenters, patientWindowSummary, '-', patColor);

style_axes(ax, ageMin, ageMax);
xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('%s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, { ...
    sprintf('6B: %s %s MFDB %s Power vs Age at %s (%s)', ...
    ageWindowModeLabel, title_case(groupSummaryStatistic), bandName, upper(electrodeLabel), upper(anesthesia)); ...
    build_window_title_line(useStationaryBins, ageWindowHalfWidthYears, ageWindowWidthYears, ageWindowStepYears)}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hCtrlCI, hCtrlWindow, hPatCI, hPatWindow], ...
    {'Controls 95% CI', sprintf('Controls %s', groupSummaryStatistic), ...
     'Patients 95% CI', sprintf('Patients %s', groupSummaryStatistic)}, ...
    'Location', 'best');

%% Step 7: Save outputs
writetable(results, fullfile(outDir, sprintf('%s_results.csv', plotFileBase)));
writetable(windowSupportTable, fullfile(outDir, sprintf('%s_window_support.csv', plotFileBase)));
save(fullfile(outDir, sprintf('%s_results.mat', plotFileBase)), ...
    'results', 'controlResults', 'patientResults', ...
    'ageWindowCenters', 'ageWindowLow', 'ageWindowHigh', 'ageWindowHalfWidthYears', ...
    'ageWindowWidthYears', 'ageWindowStepYears', 'useStationaryBins', ...
    'ageWindowModeTag', 'ageWindowModeLabel', 'groupSummaryStatistic', ...
    'controlWindowSummary', 'patientWindowSummary', ...
    'controlWindowCILow', 'controlWindowCIHigh', 'controlWindowSubjectCount', ...
    'patientWindowCILow', 'patientWindowCIHigh', 'patientWindowSubjectCount', ...
    'windowSupportTable', ...
    'subjectBootPowerDB', 'subjectBootPowerLinear', ...
    'electrodeLabel', 'bandName', 'bandRange', 'ageMin', 'ageMax', ...
    'anesthesia', 'cfg', 'validationFile');
saveas(figRaw, fullfile(outDir, sprintf('%s_plot_6A_raw_points.png', plotFileBase)));
saveas(figSlidingWindow, fullfile(outDir, sprintf('%s_plot_6B_age_window_with_ci.png', plotFileBase)));

fprintf('\nSaved outputs to %s\n', outDir);

%% Local functions
function groupLabel = standardize_group_label(groupLabel)
    if strcmpi(groupLabel, 'CP')
        groupLabel = "Patient";
    elseif strcmpi(groupLabel, 'Control')
        groupLabel = "Control";
    else
        groupLabel = string(groupLabel);
    end
end

function windowSummary = compute_age_window_summary(ageVals, powerVals, windowCenters, halfWidth, statFn)
    windowSummary = nan(numel(windowCenters), 1);
    for i = 1:numel(windowCenters)
        inWindow = is_subject_in_age_window(ageVals, windowCenters(i), halfWidth);
        if any(inWindow)
            windowSummary(i) = statFn(powerVals(inWindow), 'omitnan');
        end
    end
end

function [ciLow, ciHigh, windowSubjectCount] = compute_age_window_mfdb_ci( ...
    ageVals, bootCell, windowCenters, halfWidth, nBootGroup, statFn)

    nWindows = numel(windowCenters);
    ciLow = nan(nWindows, 1);
    ciHigh = nan(nWindows, 1);
    windowSubjectCount = zeros(nWindows, 1);

    for i = 1:nWindows
        inWindow = is_subject_in_age_window(ageVals, windowCenters(i), halfWidth);
        windowSubjectCount(i) = sum(inWindow);
        if windowSubjectCount(i) == 0
            continue;
        end

        windowBootCell = bootCell(inWindow);
        bootSummary = nan(nBootGroup, 1);
        for b = 1:nBootGroup
            drawVals = nan(windowSubjectCount(i), 1);
            for s = 1:windowSubjectCount(i)
                thisBoot = windowBootCell{s};
                idx = randi(numel(thisBoot));
                drawVals(s) = thisBoot(idx);
            end
            bootSummary(b) = statFn(drawVals, 'omitnan');
        end

        ciLow(i) = prctile(bootSummary, 2.5);
        ciHigh(i) = prctile(bootSummary, 97.5);
    end
end

function inWindow = is_subject_in_age_window(ageVals, center, halfWidth)
    inWindow = ageVals >= (center - halfWidth) & ageVals <= (center + halfWidth);
end

function style_axes(ax, ageMin, ageMax)
    ax.XLim = [ageMin ageMax];
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.Box = 'off';
    grid(ax, 'off');
end

function h = plot_valid_windows(ax, xVals, yVals, lineStyle, lineColor)
    validMask = isfinite(xVals(:)) & isfinite(yVals(:));
    if any(validMask)
        h = plot(ax, xVals(validMask), yVals(validMask), lineStyle, ...
            'Color', lineColor, 'LineWidth', 2, ...
            'Marker', 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'k');
    else
        h = plot(ax, nan, nan, lineStyle, 'Color', lineColor, 'LineWidth', 2, ...
            'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'k');
    end
end

function h = plot_valid_window_ci(ax, xVals, yVals, ciLow, ciHigh, ciColor)
    validMask = isfinite(xVals(:)) & isfinite(yVals(:)) & isfinite(ciLow(:)) & isfinite(ciHigh(:));
    if any(validMask)
        errLow = yVals(validMask) - ciLow(validMask);
        errHigh = ciHigh(validMask) - yVals(validMask);
        h = errorbar(ax, xVals(validMask), yVals(validMask), errLow, errHigh, 'LineStyle', 'none', ...
            'Color', ciColor, 'LineWidth', 1.8, 'CapSize', 8);
    else
        h = errorbar(ax, nan, nan, nan, nan, 'LineStyle', 'none', ...
            'Color', ciColor, 'LineWidth', 1.8, 'CapSize', 8);
    end
end

function out = title_case(txt)
    txt = char(string(txt));
    out = [upper(txt(1)) lower(txt(2:end))];
end

function txt = build_window_title_line(useStationaryBins, halfWidth, windowWidth, stepYears)
    if useStationaryBins
        txt = sprintf('non-overlapping bins, width %.2f yr; half-width %.2f yr; step ignored', ...
            windowWidth, halfWidth);
    else
        txt = sprintf('overlapping windows, half-width %.2f yr, step %.2f yr', ...
            halfWidth, stepYears);
    end
end
