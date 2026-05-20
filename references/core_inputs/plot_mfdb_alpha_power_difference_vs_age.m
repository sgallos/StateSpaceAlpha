%% plot_mfdb_alpha_power_difference_vs_age
% Plot MFDB alpha-power difference vs age at F3 using
% validated subject-level MFDB outputs.
%
% Difference is defined as:
%   Patient group summary alpha power - Control group summary alpha power
%
% Confidence intervals are built with the same MFDB group-bootstrap logic
% used elsewhere in the project: on each iteration, draw one bootstrap
% alpha-power value per subject within the age window, compute the
% within-group summaries, then take the group difference.
%
% Stationary non-overlapping bins are the cleaner inferential mode because
% each age point uses a distinct age interval. Sliding overlapping windows
% are the smoother trajectory-visualization mode. Both modes propagate
% upstream subject-level MFDB uncertainty into the within-window difference.

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
% Change these values to regenerate Figure 4 under different assumptions.
groupSummaryStatistic = 'mean';       % 'median' or 'mean'
stationary = 1;
ageWindowHalfWidthYears = .5;       % stationary: 0.5 or 1; sliding: 1.25 or 2
ageWindowStepYears = 0.5;             % ignored when useStationaryBins is true


% Windowing mode for the age axis. Both modes propagate MFDB uncertainty
% into the within-window group summary.
%   true  = non-overlapping bins: independent pointwise CIs, smaller per-bin n
%   false = sliding overlapping windows: smoothed trajectory, pointwise CIs only
useStationaryBins = true;

% 
% groupSummaryStatistic = 'median';
% useStationaryBins = true;
% ageWindowHalfWidthYears = 0.5;
% ageWindowStepYears = 0.5;



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

outDir = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_alpha_power_difference_vs_age');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

plotFileBase = sprintf('mfdb_%s_difference_%s_%0.1f_to_%0.1f_%s_%s_%s_hw%0.2f_step%0.2f', ...
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

V = load(validationFile, 'includedAfterValidation');
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

diffColor = [0.1 0.1 0.1];
diffCIColor = [0.75 0.80 0.86];

%% Step 5: Compute observed window summaries and MFDB-based difference CIs
controlWindowSummary = compute_age_window_summary( ...
    controlResults.AgeYears, controlResults.AlphaPower_dB, ...
    ageWindowCenters, ageWindowHalfWidthYears, groupSummaryFunction);
patientWindowSummary = compute_age_window_summary( ...
    patientResults.AgeYears, patientResults.AlphaPower_dB, ...
    ageWindowCenters, ageWindowHalfWidthYears, groupSummaryFunction);
differenceWindowSummary = patientWindowSummary - controlWindowSummary;

[differenceWindowCILow, differenceWindowCIHigh, controlWindowSubjectCount, patientWindowSubjectCount] = ...
    compute_age_window_mfdb_difference_ci( ...
    controlResults.AgeYears, controlBoot, patientResults.AgeYears, patientBoot, ...
    ageWindowCenters, ageWindowHalfWidthYears, cfg.Bgroup, groupSummaryFunction);

differenceSignificanceMask = isfinite(differenceWindowCILow) & isfinite(differenceWindowCIHigh) & ...
    ((differenceWindowCILow > 0) | (differenceWindowCIHigh < 0));

%% Step 6: Plot age-window difference with MFDB-based 95% CI
figDiff = figure('Color', 'w', 'Position', [100 100 1150 520]);
ax = axes(figDiff);
hold(ax, 'on');

for x = ageMin:1:ceil(ageMax)
    xline(ax, x, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
end

hCI = plot_valid_difference_ci(ax, ageWindowCenters, differenceWindowSummary, ...
    differenceWindowCILow, differenceWindowCIHigh, diffCIColor);
hDiff = plot_valid_windows(ax, ageWindowCenters, differenceWindowSummary, '-', diffColor);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);

style_axes(ax, ageMin, ageMax);
xlabel(ax, 'Age (years)');
ylabel(ax, sprintf('Patient - Control %s Power at %s (dB)', bandName, upper(electrodeLabel)));
title(ax, { ...
    sprintf('%s %s MFDB %s Power Difference vs Age at %s (%s)', ...
    ageWindowModeLabel, title_case(groupSummaryStatistic), bandName, upper(electrodeLabel), upper(anesthesia)); ...
    build_window_title_line(useStationaryBins, ageWindowHalfWidthYears, ageWindowWidthYears, ageWindowStepYears)}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend(ax, [hCI, hDiff, hZero], ...
    {'95% CI', sprintf('Patient - Control %s difference', groupSummaryStatistic), 'Zero line'}, ...
    'Location', 'best');

%% Step 7: Save outputs
windowDifferenceSummaryTable = table( ...
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
    differenceWindowSummary(:), ...
    differenceWindowCILow(:), ...
    differenceWindowCIHigh(:), ...
    differenceSignificanceMask(:), ...
    'VariableNames', {'AgeWindowCenter','AgeWindowLow','AgeWindowHigh', ...
    'GroupSummaryStatistic','AgeWindowMode','UseStationaryBins', ...
    'AgeWindowHalfWidthYears','AgeWindowWidthYears','AgeWindowStepYears', ...
    'ControlSubjectCount','PatientSubjectCount', ...
    'ControlSummary_dB','PatientSummary_dB','DifferenceSummary_dB', ...
    'DifferenceCILow_dB','DifferenceCIHigh_dB','CIExcludesZero'});

writetable(results, fullfile(outDir, sprintf('%s_subject_results.csv', plotFileBase)));
writetable(windowDifferenceSummaryTable, fullfile(outDir, sprintf('%s_window_difference_summary.csv', plotFileBase)));
save(fullfile(outDir, sprintf('%s_results.mat', plotFileBase)), ...
    'results', 'controlResults', 'patientResults', ...
    'ageWindowCenters', 'ageWindowLow', 'ageWindowHigh', ...
    'ageWindowHalfWidthYears', 'ageWindowWidthYears', 'ageWindowStepYears', ...
    'useStationaryBins', 'ageWindowModeTag', 'ageWindowModeLabel', ...
    'groupSummaryStatistic', ...
    'controlWindowSummary', 'patientWindowSummary', 'differenceWindowSummary', ...
    'differenceWindowCILow', 'differenceWindowCIHigh', ...
    'controlWindowSubjectCount', 'patientWindowSubjectCount', ...
    'differenceSignificanceMask', 'windowDifferenceSummaryTable', ...
    'subjectBootPowerDB', ...
    'electrodeLabel', 'bandName', 'bandRange', 'ageMin', 'ageMax', ...
    'anesthesia', 'cfg', 'validationFile');
saveas(figDiff, fullfile(outDir, sprintf('%s_plot_age_window_difference_with_ci.png', plotFileBase)));

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

function [ciLow, ciHigh, controlSubjectCount, patientSubjectCount] = compute_age_window_mfdb_difference_ci( ...
    controlAges, controlBoot, patientAges, patientBoot, windowCenters, halfWidth, nBootGroup, statFn)

    nWindows = numel(windowCenters);
    ciLow = nan(nWindows, 1);
    ciHigh = nan(nWindows, 1);
    controlSubjectCount = zeros(nWindows, 1);
    patientSubjectCount = zeros(nWindows, 1);

    for i = 1:nWindows
        inCtrl = is_subject_in_age_window(controlAges, windowCenters(i), halfWidth);
        inPat = is_subject_in_age_window(patientAges, windowCenters(i), halfWidth);
        controlSubjectCount(i) = sum(inCtrl);
        patientSubjectCount(i) = sum(inPat);

        if controlSubjectCount(i) == 0 || patientSubjectCount(i) == 0
            continue;
        end

        ctrlBootCell = controlBoot(inCtrl);
        patBootCell = patientBoot(inPat);
        bootDiff = nan(nBootGroup, 1);

        for b = 1:nBootGroup
            ctrlVals = nan(controlSubjectCount(i), 1);
            patVals = nan(patientSubjectCount(i), 1);

            for s = 1:controlSubjectCount(i)
                thisBoot = ctrlBootCell{s};
                ctrlVals(s) = thisBoot(randi(numel(thisBoot)));
            end

            for s = 1:patientSubjectCount(i)
                thisBoot = patBootCell{s};
                patVals(s) = thisBoot(randi(numel(thisBoot)));
            end

            bootDiff(b) = statFn(patVals, 'omitnan') - statFn(ctrlVals, 'omitnan');
        end

        ciLow(i) = prctile(bootDiff, 2.5);
        ciHigh(i) = prctile(bootDiff, 97.5);
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
            'Color', lineColor, 'LineWidth', 2.2, ...
            'Marker', 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'k');
    else
        h = plot(ax, nan, nan, lineStyle, 'Color', lineColor, 'LineWidth', 2.2, ...
            'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'k');
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

function h = plot_valid_difference_ci(ax, xVals, yVals, ciLow, ciHigh, ciColor)
    xVals = xVals(:);
    yVals = yVals(:);
    ciLow = ciLow(:);
    ciHigh = ciHigh(:);

    validMask = isfinite(xVals) & isfinite(yVals) & isfinite(ciLow) & isfinite(ciHigh);
    if any(validMask)
        errLow = yVals(validMask) - ciLow(validMask);
        errHigh = ciHigh(validMask) - yVals(validMask);
        h = errorbar(ax, xVals(validMask), yVals(validMask), errLow, errHigh, ...
            'LineStyle', 'none', 'Color', ciColor, 'LineWidth', 1.8, 'CapSize', 8);
    else
        h = errorbar(ax, nan, nan, nan, nan, ...
            'LineStyle', 'none', 'Color', ciColor, 'LineWidth', 1.8, 'CapSize', 8);
    end
end
