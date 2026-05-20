%% run_subject_mfdb_bootstrap_subepochs_batch
% Batch wrapper for the 4 x 30 s sub-epoch MFDB pipeline.

clear; clc;

cfg = mfdb_config();
addpath(cfg.scriptRoot);
addpath(cfg.pedsRoot);
addpath(fullfile(cfg.pedsRoot, 'src'));
addpath(genpath(fullfile(cfg.pedsRoot, 'src')));

if exist('setup_paths', 'file') == 2
    setup_paths();
end

skipExisting = false;

if ~isfile(cfg.validationFile)
    error(['Validation file not found:\n%s\n\n' ...
        'Run validate_group_mfdb_inputs.m first so this batch uses the same included cohort.'], ...
        cfg.validationFile);
end

V = load(cfg.validationFile, 'includedAfterValidation');
subjectTable = V.includedAfterValidation;
subjectTable = sortrows(subjectTable, {'groupLabel', 'ageYears', 'subjectID'});

batchOutDir = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_subepoch_batch_runs');
if ~exist(batchOutDir, 'dir')
    mkdir(batchOutDir);
end

summaryCsv = fullfile(batchOutDir, 'mfdb_subepoch_batch_summary.csv');
summaryMat = fullfile(batchOutDir, 'mfdb_subepoch_batch_summary.mat');

fprintf('Subjects selected for sub-epoch batch run: %d\n', height(subjectTable));
disp(subjectTable(:, {'subjectID', 'groupLabel', 'ageYears'}));

summary = table();
summary.subjectID = strings(0, 1);
summary.groupLabel = strings(0, 1);
summary.ageYears = zeros(0, 1);
summary.status = strings(0, 1);
summary.usedElectrodeLabel = strings(0, 1);
summary.usedFallbackElectrode = false(0, 1);
summary.nSubEpochs = zeros(0, 1);
summary.outputFile = strings(0, 1);
summary.errorMessage = strings(0, 1);

for i = 1:height(subjectTable)
    subjectID = char(subjectTable.subjectID(i));
    expectedOutDir = fullfile(cfg.subEpochOutputDir, subjectID);

    if skipExisting && isfolder(expectedOutDir)
        existing = dir(fullfile(expectedOutDir, sprintf('%s_*_mfdb_subepochs.mat', subjectID)));
        if ~isempty(existing)
            outFile = fullfile(existing(1).folder, existing(1).name);
            loaded = load(outFile, 'usedElectrodeLabel', 'usedFallbackElectrode', 'nSubEpochs');
            summary = [summary; make_summary_row(subjectTable(i, :), "skipped_existing", ...
                string(loaded.usedElectrodeLabel), logical(loaded.usedFallbackElectrode), ...
                double(loaded.nSubEpochs), string(outFile), "")]; %#ok<AGROW>
            fprintf('[%d/%d] %s skipped (existing output)\n', i, height(subjectTable), subjectID);
            continue;
        end
    end

    try
        fprintf('\n[%d/%d] Running sub-epoch MFDB for %s\n', i, height(subjectTable), subjectID);
        f3QualityPass = true;
        if ismember('f3QualityPass', subjectTable.Properties.VariableNames)
            f3QualityPass = logical(subjectTable.f3QualityPass(i));
        end

        result = run_subject_mfdb_bootstrap_subepochs_for_id(subjectID, f3QualityPass, cfg);
        summary = [summary; make_summary_row(subjectTable(i, :), "success", ...
            result.usedElectrodeLabel, result.usedFallbackElectrode, ...
            result.nSubEpochs, result.outFile, "")]; %#ok<AGROW>
    catch ME
        summary = [summary; make_summary_row(subjectTable(i, :), "error", ...
            "", false, 0, "", string(ME.message))]; %#ok<AGROW>
        fprintf('Error for %s: %s\n', subjectID, ME.message);
    end
end

writetable(summary, summaryCsv);
save(summaryMat, 'summary', 'subjectTable', 'cfg', 'skipExisting');

fprintf('\nSub-epoch batch complete.\n');
disp(summary(:, {'subjectID', 'status', 'usedElectrodeLabel', 'usedFallbackElectrode', 'nSubEpochs'}));
fprintf('Saved sub-epoch batch summary to %s\n', summaryCsv);

function row = make_summary_row(subjectRow, status, usedElectrode, usedFallback, nSubEpochs, outFile, errorMessage)
    row = table( ...
        string(subjectRow.subjectID), ...
        string(subjectRow.groupLabel), ...
        double(subjectRow.ageYears), ...
        string(status), ...
        string(usedElectrode), ...
        logical(usedFallback), ...
        double(nSubEpochs), ...
        string(outFile), ...
        string(errorMessage), ...
        'VariableNames', {'subjectID', 'groupLabel', 'ageYears', 'status', ...
        'usedElectrodeLabel', 'usedFallbackElectrode', 'nSubEpochs', ...
        'outputFile', 'errorMessage'});
end
