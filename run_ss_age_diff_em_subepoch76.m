%% run_ss_age_diff_em_subepoch76
% Diagnostic Step 4 run using all 76 sub-epoch rows directly.
%
% This is not the default manuscript path. It tests whether repeated
% within-subject observations improve identifiability of the additive
% biological variance term. Within each subject, the four 30 s sub-epochs
% are treated as repeated observations at the same age:
%
%   same subject transition:      A = I, Q = 0
%   next subject transition:      4-D IWP transition across the age gap
%
% Because same-subject transitions have zero process covariance, this
% diagnostic uses the approximate M-step for q_f0 and q_delta.

clear; clc;

%% Step 1: Locate the 76-row sub-epoch alpha table
repoRoot = fileparts(mfilename('fullpath'));
outDir = fullfile(repoRoot, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

localSubepochTableFile = fullfile(outDir, 'alpha_table_for_ssm_subepochs.csv');
upstreamSubepochTableFile = fullfile('/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS', ...
    'outputs', 'alpha_table_for_ssm_subepochs.csv');

if ~isfile(localSubepochTableFile)
    if isfile(upstreamSubepochTableFile)
        copyfile(upstreamSubepochTableFile, localSubepochTableFile);
        fprintf('Copied sub-epoch table into StateSpaceAlpha outputs:\n  %s\n', localSubepochTableFile);
    else
        error(['Sub-epoch alpha table not found.\nExpected either:\n  %s\nor:\n  %s\n\n' ...
            'Run smoke_test_alpha_extraction_subepochs from MATLAB_Multitaper_Hz_Domain_BTS first.'], ...
            localSubepochTableFile, upstreamSubepochTableFile);
    end
end

alphaTable = readtable(localSubepochTableFile, 'TextType', 'string');
fprintf('Loaded 76-row sub-epoch table: %s\n', localSubepochTableFile);

%% Step 2: Validate table shape before fitting
validate_subepoch_table(alphaTable);

fprintf('Rows: %d | Unique subjects: %d\n', height(alphaTable), numel(unique(alphaTable.subjectID)));
fprintf('CP rows: %d | Control rows: %d\n', ...
    sum(strcmpi(alphaTable.groupLabel, "CP")), sum(strcmpi(alphaTable.groupLabel, "Control")));
disp(head(sortrows(alphaTable, {'ageYears','subjectID','subepochIdx'}), 12));

%% Step 3: Run repeated-row EM diagnostic
emOptions = struct();
emOptions.maxIter = 100;
emOptions.tolerance = 1e-5;
emOptions.qInitMultipliers = [0.01 0.1 1 10 100];
emOptions.useApproximateMstep = true;

emResults = SS_age_diff_em(alphaTable, ...
    'maxIter', emOptions.maxIter, ...
    'tolerance', emOptions.tolerance, ...
    'qInitMultipliers', emOptions.qInitMultipliers, ...
    'useApproximateMstep', emOptions.useApproximateMstep, ...
    'allowRepeatedSubjectRows', true, ...
    'verbose', true);

bestResult = emResults.bestStart;
rowTrajectory = emResults.trajectory;
subjectTrajectory = keep_last_subepoch_per_subject(rowTrajectory);
multiStartSummary = build_multistart_summary(emResults);

%% Step 4: Save diagnostic outputs
trajectoryCsvFile = fullfile(outDir, 'ssm_step4_em_age_difference_trajectory_subepoch76.csv');
subjectTrajectoryCsvFile = fullfile(outDir, 'ssm_step4_em_age_difference_subject_trajectory_subepoch76.csv');
summaryCsvFile = fullfile(outDir, 'ssm_step4_em_multistart_summary_subepoch76.csv');
matFile = fullfile(outDir, 'ssm_step4_em_age_difference_results_subepoch76.mat');

writetable(rowTrajectory, trajectoryCsvFile);
writetable(subjectTrajectory, subjectTrajectoryCsvFile);
writetable(multiStartSummary, summaryCsvFile);
save(matFile, 'alphaTable', 'emResults', 'bestResult', 'rowTrajectory', ...
    'subjectTrajectory', 'multiStartSummary', 'emOptions', 'localSubepochTableFile');

fprintf('\nSaved row-level trajectory to %s\n', trajectoryCsvFile);
fprintf('Saved subject-level trajectory to %s\n', subjectTrajectoryCsvFile);
fprintf('Saved multi-start summary to %s\n', summaryCsvFile);
fprintf('Saved MAT results to %s\n', matFile);

%% Step 5: Primary diagnostic delta plot
figDelta = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figDelta);
hold(ax, 'on');

hBand = plot_ci_band(ax, subjectTrajectory.ageYears, ...
    subjectTrajectory.deltaCILow_dB, subjectTrajectory.deltaCIHigh_dB, ...
    [0.78 0.80 0.86], 0.42);
hZero = yline(ax, 0, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.4);
hDelta = plot(ax, subjectTrajectory.ageYears, subjectTrajectory.deltaMean_dB, ...
    '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 2.8);

xlabel(ax, 'Age (years)');
ylabel(ax, 'CP - Control alpha power difference (dB)');
title(ax, sprintf(['Sub-epoch 76-row diagnostic: CP-Control alpha-power difference ' ...
    '(qF0 %.3g, qDelta %.3g, sigmaBio %.3g)'], ...
    emResults.q_f0_em, emResults.q_delta_em, emResults.sigmaBio_em));
legend(ax, [hBand, hDelta, hZero], {'Delta 95% CI', 'Delta mean', 'Zero line'}, ...
    'Location', 'best');
style_ssm_axes(ax);

deltaPlotFile = fullfile(outDir, 'ssm_step4_em_delta_primary_subepoch76.png');
saveas(figDelta, deltaPlotFile);
fprintf('Saved delta plot to %s\n', deltaPlotFile);

%% Step 6: Log-likelihood diagnostic
figLogLik = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figLogLik);
hold(ax, 'on');

for i = 1:numel(emResults.startResults)
    thisStart = emResults.startResults{i};
    plot(ax, 1:numel(thisStart.logLikHistory), thisStart.logLikHistory, ...
        '-o', 'LineWidth', 1.6, 'MarkerSize', 4, ...
        'DisplayName', sprintf('start %.2g', emResults.qInitMultipliers(i)));
end

xlabel(ax, 'EM iteration');
ylabel(ax, 'Log-likelihood');
title(ax, 'Sub-epoch 76-row diagnostic: EM convergence');
legend(ax, 'Location', 'best');
style_ssm_axes(ax);

logLikPlotFile = fullfile(outDir, 'ssm_step4_em_loglik_history_subepoch76.png');
saveas(figLogLik, logLikPlotFile);
fprintf('Saved log-likelihood plot to %s\n', logLikPlotFile);

%% Step 7: Hyperparameter history diagnostic
figQ = figure('Color', 'w', 'Position', [100 100 1050 560]);
ax = axes(figQ);
hold(ax, 'on');

for i = 1:numel(emResults.startResults)
    thisStart = emResults.startResults{i};
    iterVals = 1:size(thisStart.qHistory, 1);
    semilogy(ax, iterVals, thisStart.qHistory(:, 1), '-', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q f0 start %.2g', emResults.qInitMultipliers(i)));
    semilogy(ax, iterVals, thisStart.qHistory(:, 2), '--', ...
        'LineWidth', 1.7, 'DisplayName', sprintf('q delta start %.2g', emResults.qInitMultipliers(i)));
    semilogy(ax, iterVals, thisStart.sigmaBioHistory, ':', ...
        'LineWidth', 2.0, 'DisplayName', sprintf('sigma bio start %.2g', emResults.qInitMultipliers(i)));
end

xlabel(ax, 'EM iteration');
ylabel(ax, 'Hyperparameter value');
title(ax, 'Sub-epoch 76-row diagnostic: hyperparameter histories');
legend(ax, 'Location', 'eastoutside');
set(ax, 'YScale', 'log');
style_ssm_axes(ax);

qPlotFile = fullfile(outDir, 'ssm_step4_em_hyperparameter_history_subepoch76.png');
saveas(figQ, qPlotFile);
fprintf('Saved hyperparameter plot to %s\n', qPlotFile);

%% Step 8: Console summary
fprintf('\nSub-epoch 76-row diagnostic summary\n');
fprintf('  q_f0_em: %.6g\n', emResults.q_f0_em);
fprintf('  q_delta_em: %.6g\n', emResults.q_delta_em);
fprintf('  sigmaBio_em: %.6g dB^2\n', emResults.sigmaBio_em);
fprintf('  final log-likelihood: %.6f\n', emResults.finalLogLikelihood);
fprintf('  converged: %d\n', emResults.converged);
fprintf('  best start log-likelihood monotone: %d\n', bestResult.logLikelihoodIsMonotone);

%% Local helpers
function validate_subepoch_table(T)
    requiredNames = ["subjectID","subepochIdx","groupLabel","ageYears","alpha_dB","mfdb_var"];
    names = string(T.Properties.VariableNames);
    missingNames = setdiff(requiredNames, names);
    if ~isempty(missingNames)
        error('Sub-epoch table is missing required columns: %s', strjoin(missingNames, ', '));
    end

    [G, subjectIDs] = findgroups(T.subjectID);
    rowsPerSubject = splitapply(@numel, T.subjectID, G);
    if any(rowsPerSubject ~= 4)
        badSubjects = subjectIDs(rowsPerSubject ~= 4);
        error('Expected exactly 4 sub-epochs per subject. Bad subjects: %s', strjoin(badSubjects, ', '));
    end

    for i = 1:numel(subjectIDs)
        rows = find(G == i);
        thisTable = sortrows(T(rows, :), 'subepochIdx');
        if ~isequal(thisTable.subepochIdx(:), (1:4).')
            error('Subject %s does not have subepochIdx 1:4.', subjectIDs(i));
        end
        if numel(unique(thisTable.ageYears)) ~= 1
            error('Subject %s has multiple ages.', subjectIDs(i));
        end
        if numel(unique(thisTable.groupLabel)) ~= 1
            error('Subject %s has multiple group labels.', subjectIDs(i));
        end
    end
end

function subjectTrajectory = keep_last_subepoch_per_subject(rowTrajectory)
    [G, subjectIDs] = findgroups(rowTrajectory.subjectID);
    keepMask = false(height(rowTrajectory), 1);
    for i = 1:numel(subjectIDs)
        rows = find(G == i);
        [~, orderIdx] = sort(rowTrajectory.subepochIdx(rows));
        sortedRows = rows(orderIdx);
        keepMask(sortedRows(end)) = true;
    end
    subjectTrajectory = rowTrajectory(keepMask, :);
    subjectTrajectory = sortrows(subjectTrajectory, {'ageYears','subjectID'});
end

function summaryTable = build_multistart_summary(emResults)
    nStarts = numel(emResults.startResults);
    startMultiplier = emResults.qInitMultipliers(:);
    initialQF0 = nan(nStarts, 1);
    initialQDelta = nan(nStarts, 1);
    initialSigmaBio = nan(nStarts, 1);
    finalQF0 = nan(nStarts, 1);
    finalQDelta = nan(nStarts, 1);
    finalSigmaBio = nan(nStarts, 1);
    finalLogLikelihood = nan(nStarts, 1);
    emIterations = nan(nStarts, 1);
    converged = false(nStarts, 1);
    logLikelihoodIsMonotone = false(nStarts, 1);
    hitQFloor = false(nStarts, 1);
    hitQCeiling = false(nStarts, 1);
    hitSigmaBioFloor = false(nStarts, 1);
    hitSigmaBioCeiling = false(nStarts, 1);
    isBestStart = false(nStarts, 1);

    for i = 1:nStarts
        thisStart = emResults.startResults{i};
        initialQF0(i) = thisStart.initialQF0;
        initialQDelta(i) = thisStart.initialQDelta;
        initialSigmaBio(i) = thisStart.initialSigmaBio;
        finalQF0(i) = thisStart.q_f0_em;
        finalQDelta(i) = thisStart.q_delta_em;
        finalSigmaBio(i) = thisStart.sigmaBio_em;
        finalLogLikelihood(i) = thisStart.finalLogLikelihood;
        emIterations(i) = thisStart.emIterations;
        converged(i) = thisStart.converged;
        logLikelihoodIsMonotone(i) = thisStart.logLikelihoodIsMonotone;
        hitQFloor(i) = thisStart.hitQFloor;
        hitQCeiling(i) = thisStart.hitQCeiling;
        hitSigmaBioFloor(i) = thisStart.hitSigmaBioFloor;
        hitSigmaBioCeiling(i) = thisStart.hitSigmaBioCeiling;
    end

    isBestStart(emResults.bestStartIndex) = true;
    summaryTable = table(startMultiplier, initialQF0, initialQDelta, initialSigmaBio, ...
        finalQF0, finalQDelta, finalSigmaBio, finalLogLikelihood, emIterations, ...
        converged, logLikelihoodIsMonotone, hitQFloor, hitQCeiling, ...
        hitSigmaBioFloor, hitSigmaBioCeiling, isBestStart);
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
