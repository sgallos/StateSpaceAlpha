function T = smoke_test_alpha_extraction_subepochs()
%SMOKE_TEST_ALPHA_EXTRACTION_SUBEPOCHS Build the 76-row alpha table.
%
% Each row is one subject sub-epoch. Alpha power is computed with the same
% rule as the current MFDB alpha pipeline: integrate 8-13 Hz in linear
% power, then convert the integrated value back to dB.

    cfg = mfdb_config();
    bandRange = [8 13];

    if ~isfile(cfg.validationFile)
        error(['Validation file not found:\n%s\n\n' ...
            'Run validate_group_mfdb_inputs.m first.'], cfg.validationFile);
    end

    V = load(cfg.validationFile, 'includedAfterValidation');
    registry = V.includedAfterValidation;

    nSubjects = height(registry);
    nTotal = nSubjects * cfg.subEpochCount;

    T = table('Size', [nTotal, 7], ...
        'VariableTypes', {'string', 'double', 'string', 'double', ...
        'double', 'double', 'double'}, ...
        'VariableNames', {'subjectID', 'subepochIdx', 'groupLabel', ...
        'ageYears', 'alpha_dB', 'mfdb_var', 'nFreqBins'});

    rowIdx = 0;

    for i = 1:nSubjects
        subjectID = string(registry.subjectID(i));
        subepochFile = find_subepoch_file(cfg, subjectID);

        if subepochFile == ""
            warning('Missing sub-epoch file for %s. Skipping.', subjectID);
            continue;
        end

        S = load(subepochFile, 'subEpochResults', 'groupLabel', 'ageYears');

        for subepochIdx = 1:numel(S.subEpochResults)
            thisSubEpoch = S.subEpochResults{subepochIdx};
            freq = thisSubEpoch.freq(:).';
            bandMask = freq >= bandRange(1) & freq <= bandRange(2);

            if ~any(bandMask)
                error('No alpha bins found for %s sub-epoch %d.', subjectID, subepochIdx);
            end

            originalLinear = trapz(freq(bandMask), ...
                10.^(thisSubEpoch.S_subject_original_dB(bandMask) / 10));
            alphaDB = 10 * log10(max(originalLinear, eps));

            bootLinear = trapz(freq(bandMask), ...
                10.^(thisSubEpoch.S_subject_boot_dB(:, bandMask) / 10), 2);
            bootDB = 10 * log10(max(bootLinear, eps));

            rowIdx = rowIdx + 1;
            T.subjectID(rowIdx) = subjectID;
            T.subepochIdx(rowIdx) = subepochIdx;
            T.groupLabel(rowIdx) = string(S.groupLabel);
            T.ageYears(rowIdx) = double(S.ageYears);
            T.alpha_dB(rowIdx) = alphaDB;
            T.mfdb_var(rowIdx) = var(bootDB);
            T.nFreqBins(rowIdx) = sum(bandMask);
        end
    end

    T = T(1:rowIdx, :);

    outDir = fullfile(cfg.scriptRoot, 'outputs');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    outCsv = fullfile(outDir, 'alpha_table_for_ssm_subepochs.csv');
    outMat = fullfile(outDir, 'alpha_table_for_ssm_subepochs.mat');
    writetable(T, outCsv);
    save(outMat, 'T', 'cfg', 'bandRange');

    fprintf('Wrote %d rows to %s\n', height(T), outCsv);
    fprintf('Saved MAT table to %s\n', outMat);
    disp(T(1:min(12, height(T)), :));
end

function subepochFile = find_subepoch_file(cfg, subjectID)
    subjectDir = fullfile(cfg.subEpochOutputDir, char(subjectID));
    files = dir(fullfile(subjectDir, sprintf('%s_*_mfdb_subepochs.mat', char(subjectID))));

    if isempty(files)
        subepochFile = "";
        return;
    end

    subepochFile = string(fullfile(files(1).folder, files(1).name));
end
