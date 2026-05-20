%% smoke_test_alpha_extraction
% Confirm that the new StateSpaceAlpha extraction logic reproduces the
% per-subject MFDB alpha-power values from the existing pipeline.
%
% This script does not fit OLS and does not fit the SSM. It only checks the
% input extraction step that the SSM will use.

clear; clc;

%% Step 1: Locate inputs and outputs
repoRoot = fileparts(fileparts(mfilename('fullpath')));
registryFile = fullfile(repoRoot, 'references', 'core_inputs', 'mfdb_subject_registry.csv');
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~isfile(registryFile)
    error('Registry file not found: %s', registryFile);
end

% Optional comparison target from the existing MFDB alpha-power pipeline.
% The smoke test still runs if this file is absent.
referenceAlphaCsv = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'mfdb_alpha_power_vs_age', 'mfdb_alpha_f3_10.0_to_15.5_sevo_results.csv');

bandRangeHz = [8 13];

%% Step 2: Load registry and keep included subjects only
registry = readtable(registryFile, 'TextType', 'string');

if ismember('excluded', registry.Properties.VariableNames)
    registry = registry(registry.excluded == 0, :);
end

registry = sortrows(registry, {'groupLabel', 'ageYears', 'subjectID'});
nSubjects = height(registry);

fprintf('Included subjects in registry: %d\n', nSubjects);
disp(registry(:, {'subjectID', 'groupLabel', 'ageYears', 'usedElectrodeLabel'}));

%% Step 3: Extract alpha power and MFDB variance
T = table('Size', [nSubjects, 11], ...
    'VariableTypes', {'string','string','double','double','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var', ...
    'nBoot','alpha_dB_from_dB_fields','selfCheckAbsDiff_dB', ...
    'reference_alpha_dB','referenceAbsDiff_dB','sourceFile'});

for i = 1:nSubjects
    matPath = char(registry.filePath(i));
    if ~isfile(matPath)
        error('MFDB file not found for %s: %s', registry.subjectID(i), matPath);
    end

    S = load(matPath, 'freq', ...
        'S_subject_original_linear', 'S_subject_boot_linear', ...
        'S_subject_original_dB', 'S_subject_boot_dB');

    freq = S.freq(:);
    bandMask = freq >= bandRangeHz(1) & freq <= bandRangeHz(2);
    if ~any(bandMask)
        error('No alpha-band frequency bins found for %s.', registry.subjectID(i));
    end

    % Existing pipeline convention:
    % integrate linear power over 8-13 Hz, then convert that integrated
    % alpha-band value to dB.
    alphaLinear = trapz(freq(bandMask), S.S_subject_original_linear(bandMask));
    alphaDB = 10 * log10(max(alphaLinear, eps));

    % Equivalent check using stored dB spectra: convert back to linear,
    % integrate, then convert the integrated value to dB.
    alphaLinearFromDB = trapz(freq(bandMask), 10.^(S.S_subject_original_dB(bandMask) / 10));
    alphaDBFromDBFields = 10 * log10(max(alphaLinearFromDB, eps));

    bootLinear = trapz(freq(bandMask), S.S_subject_boot_linear(:, bandMask), 2);
    bootDB = 10 * log10(max(bootLinear, eps));
    bootDB = bootDB(isfinite(bootDB));

    T.subjectID(i) = string(registry.subjectID(i));
    T.groupLabel(i) = string(registry.groupLabel(i));
    T.ageYears(i) = double(registry.ageYears(i));
    T.alpha_dB(i) = alphaDB;
    T.mfdb_var(i) = var(bootDB, 0);
    T.nBoot(i) = numel(bootDB);
    T.alpha_dB_from_dB_fields(i) = alphaDBFromDBFields;
    T.selfCheckAbsDiff_dB(i) = abs(alphaDB - alphaDBFromDBFields);
    T.reference_alpha_dB(i) = NaN;
    T.referenceAbsDiff_dB(i) = NaN;
    T.sourceFile(i) = string(matPath);
end

%% Step 4: Optional comparison against existing alpha-power CSV
if isfile(referenceAlphaCsv)
    referenceTable = readtable(referenceAlphaCsv, 'TextType', 'string');
    if all(ismember({'SubjectID', 'AlphaPower_dB'}, referenceTable.Properties.VariableNames))
        for i = 1:height(T)
            matchIdx = find(strcmpi(referenceTable.SubjectID, T.subjectID(i)), 1, 'first');
            if ~isempty(matchIdx)
                T.reference_alpha_dB(i) = referenceTable.AlphaPower_dB(matchIdx);
                T.referenceAbsDiff_dB(i) = abs(T.alpha_dB(i) - T.reference_alpha_dB(i));
            end
        end
    else
        warning('Reference CSV found, but it does not contain SubjectID and AlphaPower_dB.');
    end
else
    fprintf('Optional reference CSV not found. Skipping old-pipeline CSV comparison.\n');
end

%% Step 5: Print and save the smoke-test table
disp(T(:, {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var','nBoot', ...
    'selfCheckAbsDiff_dB','referenceAbsDiff_dB'}));

fprintf('\nAlpha dB range: %.4f to %.4f dB\n', min(T.alpha_dB), max(T.alpha_dB));
fprintf('MFDB variance range: %.6f to %.6f dB^2\n', min(T.mfdb_var), max(T.mfdb_var));
fprintf('Max self-check difference: %.3g dB\n', max(T.selfCheckAbsDiff_dB));

if any(isfinite(T.referenceAbsDiff_dB))
    fprintf('Max reference CSV difference: %.3g dB\n', max(T.referenceAbsDiff_dB, [], 'omitnan'));
end

outCsv = fullfile(outDir, 'alpha_table_for_ssm.csv');
writetable(T, outCsv);
fprintf('Saved alpha table for SSM to %s\n', outCsv);

%% Step 6: Figure 1, smoke-test scatter
fig = figure('Color', 'w', 'Position', [100 100 850 500]);
ax = axes(fig);
hold(ax, 'on');

controlMask = strcmpi(T.groupLabel, 'Control');
cpMask = strcmpi(T.groupLabel, 'CP') | strcmpi(T.groupLabel, 'Patient');

hControl = scatter(ax, T.ageYears(controlMask), T.alpha_dB(controlMask), 70, ...
    'MarkerFaceColor', [0.9 0.45 0.1], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
hCP = scatter(ax, T.ageYears(cpMask), T.alpha_dB(cpMask), 70, ...
    'MarkerFaceColor', [0.2 0.45 0.85], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);

xlabel(ax, 'Age (yr)');
ylabel(ax, 'Alpha power (dB)');
title(ax, 'Figure 1: smoke-test extracted MFDB alpha power vs age');
legend(ax, [hControl, hCP], {'Control', 'CP'}, 'Location', 'best');
ax.FontName = 'Helvetica';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.TickDir = 'out';
box(ax, 'off');

plotFile = fullfile(outDir, 'smoke_test_alpha_vs_age.png');
saveas(fig, plotFile);
fprintf('Saved smoke-test plot to %s\n', plotFile);

%% Step 7: Figure 2, comparison against existing pipeline scatter
hasReference = any(isfinite(T.reference_alpha_dB));

if hasReference
    figCompare = figure('Color', 'w', 'Position', [100 100 1150 500]);
    tiledlayout(figCompare, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Left panel: same alpha-vs-age scatter, with existing pipeline values
    % plotted as open markers and smoke-test values plotted as filled markers.
    axOverlay = nexttile;
    hold(axOverlay, 'on');

    hRefControl = scatter(axOverlay, T.ageYears(controlMask), T.reference_alpha_dB(controlMask), 110, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.9 0.45 0.1], 'LineWidth', 1.8);
    hRefCP = scatter(axOverlay, T.ageYears(cpMask), T.reference_alpha_dB(cpMask), 110, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.2 0.45 0.85], 'LineWidth', 1.8);
    hSmokeControl = scatter(axOverlay, T.ageYears(controlMask), T.alpha_dB(controlMask), 45, ...
        'MarkerFaceColor', [0.9 0.45 0.1], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    hSmokeCP = scatter(axOverlay, T.ageYears(cpMask), T.alpha_dB(cpMask), 45, ...
        'MarkerFaceColor', [0.2 0.45 0.85], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

    xlabel(axOverlay, 'Age (yr)');
    ylabel(axOverlay, 'Alpha power (dB)');
    title(axOverlay, 'Alpha-vs-age overlay');
    legend(axOverlay, [hRefControl, hSmokeControl, hRefCP, hSmokeCP], ...
        {'Existing Control', 'Smoke Control', 'Existing CP', 'Smoke CP'}, ...
        'Location', 'best');
    style_smoke_axes(axOverlay);

    % Right panel: numerical agreement with the old pipeline.
    axIdentity = nexttile;
    hold(axIdentity, 'on');

    allVals = [T.reference_alpha_dB(:); T.alpha_dB(:)];
    plotLimits = [floor(min(allVals, [], 'omitnan') - 1), ceil(max(allVals, [], 'omitnan') + 1)];

    plot(axIdentity, plotLimits, plotLimits, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.2);
    scatter(axIdentity, T.reference_alpha_dB(controlMask), T.alpha_dB(controlMask), 70, ...
        'MarkerFaceColor', [0.9 0.45 0.1], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    scatter(axIdentity, T.reference_alpha_dB(cpMask), T.alpha_dB(cpMask), 70, ...
        'MarkerFaceColor', [0.2 0.45 0.85], 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);

    xlabel(axIdentity, 'Existing pipeline alpha power (dB)');
    ylabel(axIdentity, 'Smoke-test alpha power (dB)');
    title(axIdentity, sprintf('Agreement check, max abs diff = %.3g dB', ...
        max(T.referenceAbsDiff_dB, [], 'omitnan')));
    xlim(axIdentity, plotLimits);
    ylim(axIdentity, plotLimits);
    axis(axIdentity, 'square');
    style_smoke_axes(axIdentity);

    comparePlotFile = fullfile(outDir, 'figure_2_smoke_test_vs_existing_pipeline.png');
    saveas(figCompare, comparePlotFile);
    fprintf('Saved Figure 2 comparison plot to %s\n', comparePlotFile);
else
    fprintf('Figure 2 comparison was skipped because no reference alpha values were available.\n');
end

function style_smoke_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
end
