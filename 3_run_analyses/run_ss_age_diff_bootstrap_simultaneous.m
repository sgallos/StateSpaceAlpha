%% run_ss_age_diff_bootstrap_simultaneous
% Step 5 entry point for MFDB bootstrap simultaneous confidence bands.
%
% Debug first with nBootstrapIterations = 50. Once the cache and figures
% look right, change nBootstrapIterations to 2000 for the meeting analysis.

clear; clc;

%% Step 1: User-facing analysis parameters
% Use alpha_table_for_ssm.csv for the current 120 s MFDB bootstrap.
% A collapsed 4 x 30 s table can be fit by Step 3/4 now, but Step 5 still
% needs a sub-epoch-aware bootstrap cache before it should be used there.
alphaTableFileName = 'alpha_table_for_ssm.csv';
nBootstrapIterations = 50;      % use 2000 for the final run
qSource = "em";                 % "manual" or "em"
bootstrapMode = "two-stage";    % "two-stage" or "mfdb-only"
manualQScaleF0 = 1;
manualQScaleDelta = 1;
manualProcessNoiseQF0 = [];
manualProcessNoiseQDelta = [];
manualBiologicalVariance = 2.451431;
bandRange = [8 13];
rngSeed = 1;
rebuildCache = false;

%% Step 2: Locate inputs
repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));
alphaTableFile = fullfile(outDir, alphaTableFileName);
registryFile = fullfile(repoRoot, 'references', 'core_inputs', 'mfdb_subject_registry.csv');
cacheFile = fullfile(outDir, 'boot_alpha_dB_cache.mat');

if ~isfile(alphaTableFile)
    error(['Alpha table not found:\n%s\n\n' ...
        'Run this first:\n' ...
        '  cd(''%s'')\n' ...
        '  smoke_test_alpha_extraction'], alphaTableFile, repoRoot);
end

alphaTable = readtable(alphaTableFile, 'TextType', 'string');
fprintf('Loaded alpha table: %s\n', alphaTableFile);

%% Step 3: Choose fixed q values for bootstrap refits
qSource = lower(qSource);
switch qSource
    case "manual"
        processNoiseQF0 = manualProcessNoiseQF0;
        processNoiseQDelta = manualProcessNoiseQDelta;
        qScaleF0 = manualQScaleF0;
        qScaleDelta = manualQScaleDelta;
        biologicalVariance = manualBiologicalVariance;
        qDescription = sprintf('manual q scales: f0 %.3g, delta %.3g', qScaleF0, qScaleDelta);

    case "em"
        step4MatFile = fullfile(outDir, 'ssm_step4_em_age_difference_results.mat');
        if ~isfile(step4MatFile)
            error('Step 4 EM results not found. Run run_ss_age_diff_em first or set qSource = "manual".');
        end

        S = load(step4MatFile, 'emResults');
        processNoiseQF0 = S.emResults.q_f0_em;
        processNoiseQDelta = S.emResults.q_delta_em;
        if isfield(S.emResults, 'sigmaBio_em')
            biologicalVariance = S.emResults.sigmaBio_em;
        else
            error('Step 4 EM results do not include sigmaBio_em. Rerun run_ss_age_diff_em.');
        end
        qScaleF0 = 1;
        qScaleDelta = 1;
        qDescription = sprintf('EM q values: f0 %.3g, delta %.3g, sigmaBio %.3g', ...
            processNoiseQF0, processNoiseQDelta, biologicalVariance);

        if isfield(S.emResults, 'hitQFloor') && ...
                (S.emResults.hitQFloor || S.emResults.hitQCeiling || ...
                S.emResults.hitSigmaBioFloor || S.emResults.hitSigmaBioCeiling)
            warning(['Step 4 EM q values hit a boundary. The bootstrap will run, but these ' ...
                'hyperparameters should be treated cautiously.']);
        end

    otherwise
        error('Unknown qSource "%s". Use "manual" or "em".', qSource);
end

fprintf('Step 5 q source: %s\n', qDescription);

%% Step 4: Run MFDB bootstrap simultaneous-band wrapper
bootResults = SS_age_diff_bootstrap_simultaneous(alphaTable, ...
    'B', nBootstrapIterations, ...
    'bandRange', bandRange, ...
    'registryFile', registryFile, ...
    'cacheFile', cacheFile, ...
    'rebuildCache', rebuildCache, ...
    'rngSeed', rngSeed, ...
    'qScaleF0', qScaleF0, ...
    'qScaleDelta', qScaleDelta, ...
    'processNoiseQF0', processNoiseQF0, ...
    'processNoiseQDelta', processNoiseQDelta, ...
    'biologicalVariance', biologicalVariance, ...
    'bootstrapMode', bootstrapMode, ...
    'verbose', true);

%% Step 5: Save Step 5 outputs
step5MatFile = fullfile(outDir, 'ssm_step5_bootstrap_simultaneous_results.mat');
step5BandCsvFile = fullfile(outDir, 'ssm_step5_bootstrap_band_summary.csv');
step5DistributionCsvFile = fullfile(outDir, 'ssm_step5_bootstrap_delta_distribution.csv');

bandSummary = bootResults.bandSummary;
deltaDistributionTable = make_delta_distribution_table(bootResults);

writetable(bandSummary, step5BandCsvFile);
writetable(deltaDistributionTable, step5DistributionCsvFile);
save(step5MatFile, 'bootResults', 'bandSummary', 'deltaDistributionTable', ...
    'nBootstrapIterations', 'qSource', 'qDescription', 'bootstrapMode', ...
    'biologicalVariance', 'bandRange', 'rngSeed', ...
    'alphaTableFileName', 'alphaTableFile');

fprintf('\nSaved Step 5 MAT results to %s\n', step5MatFile);
fprintf('Saved Step 5 band summary to %s\n', step5BandCsvFile);
fprintf('Saved Step 5 delta bootstrap distribution to %s\n', step5DistributionCsvFile);

%% Step 6: Figure 1, primary simultaneous band
figSim = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figSim);
hold(ax, 'on');

hSimBand = plot_ci_band(ax, bootResults.ageGrid, ...
    bootResults.deltaSimultaneousLow, bootResults.deltaSimultaneousHigh, ...
    [0.70 0.75 0.86], 0.42);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hOriginal = plot(ax, bootResults.ageGrid, ...
    bootResults.originalFit.trajectory.deltaMean_dB(:).', ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

highlight_significant_windows(ax, bootResults.ageGrid, bootResults.simultaneousSigMask);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, sprintf('CP-Control alpha-power difference, simultaneous 95%% band (B = %d)', ...
    bootResults.diagnostics.nValidBootstraps));
legend(ax, [hSimBand, hOriginal, hZero], ...
    {'Simultaneous 95% band', 'Original fitted delta', 'Zero line'}, 'Location', 'best');
style_ssm_axes(ax);

simPlotFile = fullfile(outDir, 'ssm_step5_delta_simultaneous_band.png');
saveas(figSim, simPlotFile);
fprintf('Saved Step 5 simultaneous-band plot to %s\n', simPlotFile);

%% Step 7: Figure 2, pointwise vs simultaneous bands
figCompare = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figCompare);
hold(ax, 'on');

hSimBand = plot_ci_band(ax, bootResults.ageGrid, ...
    bootResults.deltaSimultaneousLow, bootResults.deltaSimultaneousHigh, ...
    [0.55 0.62 0.78], 0.30);
hPointBand = plot_ci_band(ax, bootResults.ageGrid, ...
    bootResults.deltaPointwiseLow, bootResults.deltaPointwiseHigh, ...
    [0.82 0.84 0.90], 0.48);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hOriginal = plot(ax, bootResults.ageGrid, ...
    bootResults.originalFit.trajectory.deltaMean_dB(:).', ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, sprintf('Pointwise vs simultaneous bands (T* = %.2f)', bootResults.diagnostics.tStar));
legend(ax, [hSimBand, hPointBand, hOriginal, hZero], ...
    {'Simultaneous 95% band', 'Pointwise 95% band', 'Original fitted delta', 'Zero line'}, ...
    'Location', 'best');
style_ssm_axes(ax);

comparePlotFile = fullfile(outDir, 'ssm_step5_pointwise_vs_simultaneous_band.png');
saveas(figCompare, comparePlotFile);
fprintf('Saved Step 5 pointwise-vs-simultaneous plot to %s\n', comparePlotFile);

%% Step 8: Figure 3, bootstrap trajectory diagnostic
figDiag = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figDiag);
hold(ax, 'on');

validRows = find(bootResults.validBootstrapMask);
nDiagnosticCurves = min(20, numel(validRows));
diagnosticRows = validRows(randperm(numel(validRows), nDiagnosticCurves));

for i = 1:numel(diagnosticRows)
    plot(ax, bootResults.ageGrid, bootResults.deltaBootstrap(diagnosticRows(i), :), ...
        '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hOriginal = plot(ax, bootResults.ageGrid, ...
    bootResults.originalFit.trajectory.deltaMean_dB(:).', ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, 'Bootstrap trajectory diagnostic');
legend(ax, [hOriginal, hZero], {'Original fitted delta', 'Zero line'}, 'Location', 'best');
style_ssm_axes(ax);

diagnosticPlotFile = fullfile(outDir, 'ssm_step5_bootstrap_trajectory_diagnostic.png');
saveas(figDiag, diagnosticPlotFile);
fprintf('Saved Step 5 bootstrap trajectory diagnostic plot to %s\n', diagnosticPlotFile);

%% Step 9: Console summary
fprintf('\n=== STEP 5 SUMMARY ===\n');
fprintf('Valid bootstrap fits: %d / %d\n', ...
    bootResults.diagnostics.nValidBootstraps, bootResults.diagnostics.nRequestedBootstraps);
fprintf('T_star simultaneous critical value: %.3f\n', bootResults.diagnostics.tStar);
fprintf('Max abs bootstrap bias: %.3f dB\n', bootResults.diagnostics.maxAbsBootstrapBias);
fprintf('Pointwise band excludes zero anywhere: %d\n', ...
    bootResults.diagnostics.pointwiseAnyExcludesZero);
fprintf('Simultaneous band excludes zero anywhere: %d\n', ...
    bootResults.diagnostics.simultaneousAnyExcludesZero);

%% Local helpers
function T = make_delta_distribution_table(bootResults)
    deltaBootstrap = bootResults.deltaBootstrap;
    nAge = numel(bootResults.ageGrid);
    varNames = strings(1, nAge);
    for i = 1:nAge
        varNames(i) = sprintf('age_%0.3f', bootResults.ageGrid(i));
        varNames(i) = replace(varNames(i), '.', 'p');
    end
    T = array2table(deltaBootstrap, 'VariableNames', cellstr(varNames));
    T.bootstrapIteration = (1:size(deltaBootstrap, 1)).';
    T = movevars(T, 'bootstrapIteration', 'Before', 1);
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

function highlight_significant_windows(ax, ageGrid, sigMask)
    if ~any(sigMask)
        return;
    end

    yLimits = ylim(ax);
    runs = find_logical_runs(sigMask(:));
    for r = 1:size(runs, 1)
        xLow = ageGrid(runs(r, 1));
        xHigh = ageGrid(runs(r, 2));
        hSig = patch(ax, [xLow xHigh xHigh xLow], ...
            [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            [1.00 0.72 0.72], 'EdgeColor', 'none', 'FaceAlpha', 0.18, ...
            'HandleVisibility', 'off');
        uistack(hSig, 'bottom');
    end
end

function runs = find_logical_runs(mask)
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    runs = [starts ends];
end

function style_ssm_axes(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 12;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'off');
end
