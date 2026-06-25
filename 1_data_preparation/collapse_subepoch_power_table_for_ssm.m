function collapsedTable = collapse_subepoch_power_table_for_ssm()
%COLLAPSE_SUBEPOCH_POWER_TABLE_FOR_SSM Collapse 4 x 30 s power rows to subjects.
%
% Input:
%   outputs/power_table_for_ssm_subepochs.csv       (76 rows)
%
% Output:
%   outputs/power_table_for_ssm_subepochs_collapsed.csv  (19 rows)
%   outputs/power_table_for_ssm_subepochs_collapsed.mat
%
% The primary scalar-r SSM uses the *_dB columns and estimates one total
% observation-noise variance r by EM. The variance/SD columns saved here are
% diagnostics and are not used by that scalar-r workflow.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');

inputCsv = fullfile(outDir, 'power_table_for_ssm_subepochs.csv');
collapsedCsv = fullfile(outDir, 'power_table_for_ssm_subepochs_collapsed.csv');
collapsedMat = fullfile(outDir, 'power_table_for_ssm_subepochs_collapsed.mat');

if ~isfile(inputCsv)
    error(['Multi-band sub-epoch table not found:\n%s\n\n' ...
        'Run extract_band_power_subepochs first.'], inputCsv);
end

subepochTable = readtable(inputCsv, 'TextType', 'string');
bands = ["delta", "theta", "alpha"];
validate_input_table(subepochTable, bands);

subepochTable = sortrows(subepochTable, {'groupLabel','ageYears','subjectID','subepochIdx'});
[G, subjectIDs] = findgroups(subepochTable.subjectID);
nSubjects = numel(subjectIDs);

collapsedTable = table();
collapsedTable.subjectID = subjectIDs;
collapsedTable.groupLabel = splitapply(@(x) x(1), subepochTable.groupLabel, G);
collapsedTable.ageYears = splitapply(@(x) x(1), subepochTable.ageYears, G);
collapsedTable.nSubEpochs = splitapply(@numel, subepochTable.subepochIdx, G);

for bandIdx = 1:numel(bands)
    bandName = bands(bandIdx);
    powerColumn = bandName + "_dB";
    internalVarColumn = bandName + "_mfdb_var";

    subjectMean = splitapply(@(x) mean(x, 'omitnan'), subepochTable.(powerColumn), G);
    subjectSD = splitapply(@(x) std(x, 0, 'omitnan'), subepochTable.(powerColumn), G);
    subjectRange = splitapply(@(x) max(x) - min(x), subepochTable.(powerColumn), G);
    empiricalSemVar = splitapply(@(x) var(x, 0, 'omitnan') / numel(x), ...
        subepochTable.(powerColumn), G);
    meanInternalVar = splitapply(@(x) mean(x, 'omitnan'), ...
        subepochTable.(internalVarColumn), G);

    collapsedTable.(powerColumn) = subjectMean;
    collapsedTable.(bandName + "_empiricalSemVariance") = max(empiricalSemVar, eps);
    collapsedTable.(bandName + "_withinSubjectSD_dB") = subjectSD;
    collapsedTable.(bandName + "_withinSubjectRange_dB") = subjectRange;
    collapsedTable.(bandName + "_meanSubepochMfdbVariance") = meanInternalVar;
    collapsedTable.(bandName + "_meanSubepochMfdbSD_dB") = sqrt(meanInternalVar);
end

writetable(collapsedTable, collapsedCsv);
save(collapsedMat, 'collapsedTable', 'subepochTable', 'bands');

fprintf('Loaded %d sub-epoch rows from %s\n', height(subepochTable), inputCsv);
fprintf('Collapsed to %d subject rows.\n', height(collapsedTable));
fprintf('Saved collapsed multi-band table to %s\n', collapsedCsv);
disp(collapsedTable(:, {'subjectID','groupLabel','ageYears','delta_dB','theta_dB','alpha_dB'}));
end

function validate_input_table(subepochTable, bands)
requiredColumns = ["subjectID", "subepochIdx", "groupLabel", "ageYears"];
for bandIdx = 1:numel(bands)
    bandName = bands(bandIdx);
    requiredColumns = [requiredColumns, bandName + "_dB", bandName + "_mfdb_var"]; %#ok<AGROW>
end

missingColumns = setdiff(requiredColumns, string(subepochTable.Properties.VariableNames));
if ~isempty(missingColumns)
    error('Input table is missing required columns: %s', strjoin(missingColumns, ', '));
end
end
