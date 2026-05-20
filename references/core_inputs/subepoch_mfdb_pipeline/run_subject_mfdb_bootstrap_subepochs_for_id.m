function result = run_subject_mfdb_bootstrap_subepochs_for_id(subjectID, f3QualityPass, cfg)
%RUN_SUBJECT_MFDB_BOOTSTRAP_SUBEPOCHS_FOR_ID Run MFDB on 30 s F3 sub-epochs.
%
% This is a parallel pipeline to run_subject_mfdb_bootstrap_for_id.m. It
% keeps the original 120 s MFDB output untouched, splits the same F3 epoch
% into four non-overlapping 30 s sub-epochs, and saves all sub-epoch MFDB
% spectra in one subject-level .mat file.

    if nargin < 2 || isempty(f3QualityPass)
        f3QualityPass = true;
    end
    if nargin < 3 || isempty(cfg)
        cfg = mfdb_config();
    end

    addpath(cfg.scriptRoot);
    addpath(cfg.pedsRoot);
    addpath(fullfile(cfg.pedsRoot, 'src'));
    addpath(genpath(fullfile(cfg.pedsRoot, 'src')));

    if exist('setup_paths', 'file') == 2
        setup_paths();
    end

    if ~exist(cfg.subEpochOutputDir, 'dir')
        mkdir(cfg.subEpochOutputDir);
    end

    fprintf('Subject: %s\n', subjectID);
    fprintf('Requested electrode: %s\n', cfg.electrodeLabel);
    fprintf('Sub-epochs: %d x %.1f s | TW = %g, K = %g, Bsubj = %d, fpass = [%g %g]\n', ...
        cfg.subEpochCount, cfg.subEpochDurationSec, cfg.subEpochTW, ...
        cfg.subEpochK, cfg.Bsubj, cfg.fpass(1), cfg.fpass(2));

    if exist('DropboxUtils', 'class') == 8
        dropboxBase = DropboxUtils.getDropboxBase();
        dataDir = fullfile(dropboxBase, 'Sebastian Gallo', ...
            'Preprocessed data', 'Sevofluorane');
    else
        dataDir = get_preprocessed_data_dir();
    end

    files = find_subject_data_files(dataDir, {subjectID});
    if isempty(files)
        error('No focus file found for %s in %s', subjectID, dataDir);
    end

    dataPath = files(1).path;
    fprintf('Loading %s\n', dataPath);

    D = load(dataPath);
    assert(isfield(D, 'data') && isfield(D, 'HDR'), ...
        'Expected fields data and HDR in %s', dataPath);

    rawData = D.data;
    if size(rawData, 1) < size(rawData, 2)
        data = rawData';
    else
        data = rawData;
    end

    if size(data, 1) < size(data, 2)
        error('Data must be samples x channels after orientation fix.');
    end

    labels = string(D.HDR.label(:));
    Fs = D.HDR.frequency(1);
    N = size(data, 1);
    T = N / Fs;

    requestedElectrodeLabel = string(cfg.electrodeLabel);
    fallbackElectrodeLabel = "F4";

    requestedIdx = find(strcmpi(labels, requestedElectrodeLabel), 1, 'first');
    fallbackIdx = find(strcmpi(labels, fallbackElectrodeLabel), 1, 'first');

    usedFallbackElectrode = false;
    electrodeSelectionNote = "";

    if ~isempty(requestedIdx)
        if ~f3QualityPass
            error('F3 is present but its quality flag is false for %s. Exclude this subject.', subjectID);
        end
        usedElectrodeLabel = requestedElectrodeLabel;
        usedElectrodeIndex = requestedIdx;
        electrodeSelectionNote = "Used F3.";
    elseif ~isempty(fallbackIdx)
        usedElectrodeLabel = fallbackElectrodeLabel;
        usedElectrodeIndex = fallbackIdx;
        usedFallbackElectrode = true;
        electrodeSelectionNote = "F3 unavailable. Used F4 fallback.";
        warning('run_subject_mfdb_bootstrap_subepochs:F4Fallback', ...
            'F3 not found for %s. Using F4 instead.', subjectID);
    else
        error('Neither F3 nor F4 was found for %s.', subjectID);
    end

    signalUsedElectrode = data(:, usedElectrodeIndex);
    samplesPerSubEpoch = round(cfg.subEpochDurationSec * Fs);
    nSubEpochs = cfg.subEpochCount;
    expectedSamples = samplesPerSubEpoch * nSubEpochs;

    if numel(signalUsedElectrode) < expectedSamples
        error(['Signal too short for %d sub-epochs of %d samples each. ' ...
            'Got %d samples.'], nSubEpochs, samplesPerSubEpoch, ...
            numel(signalUsedElectrode));
    end

    if numel(signalUsedElectrode) ~= expectedSamples
        warning('Trimming %s from %d to %d samples for exactly %d x %.1f s sub-epochs.', ...
            subjectID, numel(signalUsedElectrode), expectedSamples, ...
            nSubEpochs, cfg.subEpochDurationSec);
    end
    signalForSubEpochs = signalUsedElectrode(1:expectedSamples);

    groupLabel = infer_group_label(subjectID);
    ageYears = lookup_subject_age_years(subjectID, cfg.anesthesia);

    fprintf('Using %s index %d\n', usedElectrodeLabel, usedElectrodeIndex);
    fprintf('Group: %s | AgeYears: %.4f | Fs: %.6f Hz\n', groupLabel, ageYears, Fs);

    subEpochResults = cell(nSubEpochs, 1);
    for subEpochIdx = 1:nSubEpochs
        sampleStart = (subEpochIdx - 1) * samplesPerSubEpoch + 1;
        sampleEnd = subEpochIdx * samplesPerSubEpoch;
        signalThisSubEpoch = signalForSubEpochs(sampleStart:sampleEnd);

        fprintf('  Sub-epoch %d/%d: samples %d-%d\n', ...
            subEpochIdx, nSubEpochs, sampleStart, sampleEnd);

        subEpochResults{subEpochIdx} = compute_subepoch_mfdb( ...
            signalThisSubEpoch, Fs, cfg, subEpochIdx, sampleStart, sampleEnd);
    end

    subjectOutDir = fullfile(cfg.subEpochOutputDir, subjectID);
    if ~exist(subjectOutDir, 'dir')
        mkdir(subjectOutDir);
    end

    outFile = fullfile(subjectOutDir, sprintf('%s_%s_mfdb_subepochs.mat', ...
        subjectID, upper(char(usedElectrodeLabel))));

    metadata = struct();
    metadata.subjectID = subjectID;
    metadata.groupLabel = groupLabel;
    metadata.dataPath = dataPath;
    metadata.filePath = outFile;
    metadata.requestedElectrodeLabel = requestedElectrodeLabel;
    metadata.usedElectrodeLabel = usedElectrodeLabel;
    metadata.usedElectrodeIndex = usedElectrodeIndex;
    metadata.usedFallbackElectrode = usedFallbackElectrode;
    metadata.electrodeSelectionNote = electrodeSelectionNote;
    metadata.f3QualityPass = f3QualityPass;
    metadata.ageYears = ageYears;
    metadata.Fs = Fs;
    metadata.T = T;
    metadata.N = N;
    metadata.epochDurationSec = cfg.epochDurationSec;
    metadata.subEpochDurationSec = cfg.subEpochDurationSec;
    metadata.subEpochCount = cfg.subEpochCount;
    metadata.samplesPerSubEpoch = samplesPerSubEpoch;
    metadata.subEpochTW = cfg.subEpochTW;
    metadata.subEpochK = cfg.subEpochK;
    metadata.Bsubj = cfg.Bsubj;
    metadata.Bgroup = cfg.Bgroup;
    metadata.fpass = cfg.fpass;
    metadata.method = cfg.method;

    save(outFile, ...
        'subjectID', 'groupLabel', 'ageYears', 'subEpochResults', ...
        'nSubEpochs', 'Fs', 'T', 'N', 'samplesPerSubEpoch', ...
        'requestedElectrodeLabel', 'usedElectrodeLabel', ...
        'usedElectrodeIndex', 'usedFallbackElectrode', ...
        'electrodeSelectionNote', 'metadata', '-v7.3');

    fprintf('Saved sub-epoch MFDB output to %s\n', outFile);

    result = struct();
    result.subjectID = string(subjectID);
    result.groupLabel = string(groupLabel);
    result.outFile = string(outFile);
    result.usedElectrodeLabel = string(usedElectrodeLabel);
    result.usedFallbackElectrode = logical(usedFallbackElectrode);
    result.ageYears = ageYears;
    result.Fs = Fs;
    result.T = T;
    result.N = N;
    result.nSubEpochs = nSubEpochs;
    result.samplesPerSubEpoch = samplesPerSubEpoch;
end

function subEpoch = compute_subepoch_mfdb(signalThisSubEpoch, Fs, cfg, subEpochIdx, sampleStart, sampleEnd)
    [S_mtm, F_mtm, X_mtm, ~] = mtm(signalThisSubEpoch, ...
        cfg.subEpochTW, cfg.subEpochK, Fs, cfg.method);

    nfft = size(X_mtm, 1);
    if mod(nfft, 2) ~= 0
        error('nfft must be even. Got %d.', nfft);
    end

    positiveIdx = (1:(nfft / 2))';
    freqFull = F_mtm(positiveIdx);
    freqMask = freqFull >= cfg.fpass(1) & freqFull <= cfg.fpass(2);
    freq = freqFull(freqMask);

    if isempty(freq)
        error('No frequency bins remain after fpass [%g %g].', ...
            cfg.fpass(1), cfg.fpass(2));
    end

    S_subject_original_linear = S_mtm(positiveIdx, 1);
    S_subject_original_linear = S_subject_original_linear(freqMask).';

    X_boot = complex(zeros(nfft / 2, cfg.subEpochK, cfg.Bsubj));
    scale = sqrt(S_mtm(positiveIdx, 1) / 2);
    scale(scale == 0) = eps;

    residuals = zeros(nfft, cfg.subEpochK);
    residuals(1:(nfft / 2), :) = real(X_mtm(positiveIdx, :, 1)) ./ scale;
    residuals((nfft / 2 + 1):nfft, :) = imag(X_mtm(positiveIdx, :, 1)) ./ scale;

    residualMean = mean(residuals, 1);
    residualStd = sqrt(mean((residuals - residualMean).^2, 1));
    residualStd(residualStd == 0) = eps;
    residualNorm = (residuals - residualMean) ./ residualStd;

    bootsam = randi(nfft, [nfft, cfg.Bsubj]);

    for taperIdx = 1:cfg.subEpochK
        residualBoot = residualNorm(:, taperIdx);
        residualBoot = residualBoot(bootsam);
        X_boot(:, taperIdx, :) = residualBoot(1:(nfft / 2), :) .* scale + ...
            1i * residualBoot((nfft / 2 + 1):nfft, :) .* scale;
    end

    PSD_boot = squeeze(mean(abs(X_boot).^2, 2));
    S_subject_boot_linear = PSD_boot(freqMask, :).';

    subEpoch = struct();
    subEpoch.subEpochIdx = subEpochIdx;
    subEpoch.sampleStart = sampleStart;
    subEpoch.sampleEnd = sampleEnd;
    subEpoch.timeStartSec = (sampleStart - 1) / Fs;
    subEpoch.timeEndSec = sampleEnd / Fs;
    subEpoch.durationSec = numel(signalThisSubEpoch) / Fs;
    subEpoch.N = numel(signalThisSubEpoch);
    subEpoch.Fs = Fs;
    subEpoch.TW = cfg.subEpochTW;
    subEpoch.K = cfg.subEpochK;
    subEpoch.Bsubj = cfg.Bsubj;
    subEpoch.fpass = cfg.fpass;
    subEpoch.freq = freq(:).';
    subEpoch.S_subject_original_linear = S_subject_original_linear;
    subEpoch.S_subject_boot_linear = S_subject_boot_linear;
    subEpoch.S_subject_original_dB = 10 * log10(max(S_subject_original_linear, eps));
    subEpoch.S_subject_boot_dB = 10 * log10(max(S_subject_boot_linear, eps));
end

function groupLabel = infer_group_label(subjectID)
    if startsWith(string(subjectID), "CS", 'IgnoreCase', true)
        groupLabel = 'CP';
    elseif startsWith(string(subjectID), "CN", 'IgnoreCase', true)
        groupLabel = 'Control';
    else
        error('Cannot infer group label from subjectID %s.', subjectID);
    end
end

function ageYears = lookup_subject_age_years(subjectID, anesthesia)
    meta = load_subject_metadata(anesthesia);

    controlIDs = string(meta.controls.ids(:));
    patientIDs = string(meta.patients.ids(:));
    controlAges = [meta.controls.ages{:}]';
    patientAges = [meta.patients.ages{:}]';

    idx = find(strcmpi(controlIDs, string(subjectID)), 1, 'first');
    if ~isempty(idx)
        ageYears = controlAges(idx);
        return;
    end

    idx = find(strcmpi(patientIDs, string(subjectID)), 1, 'first');
    if ~isempty(idx)
        ageYears = patientAges(idx);
        return;
    end

    ageYears = NaN;
end
