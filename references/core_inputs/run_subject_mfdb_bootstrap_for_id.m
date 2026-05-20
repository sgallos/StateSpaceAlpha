function result = run_subject_mfdb_bootstrap_for_id(subjectID, f3QualityPass, cfg)
%RUN_SUBJECT_MFDB_BOOTSTRAP_FOR_ID Run subject-level MFDB bootstrap for one subject.

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

    if ~exist(cfg.outputDir, 'dir')
        mkdir(cfg.outputDir);
    end

    fprintf('Subject: %s\n', subjectID);
    fprintf('Electrode: %s\n', cfg.electrodeLabel);
    fprintf('TW = %g, K = %g, Bsubj = %d, fpass = [%g %g]\n', ...
        cfg.TW, cfg.K, cfg.Bsubj, cfg.fpass(1), cfg.fpass(2));

    if exist('DropboxUtils', 'class') == 8
        dropboxBase = DropboxUtils.getDropboxBase();
        dataDir = fullfile(dropboxBase, 'Sebastian Gallo', 'Preprocessed data', 'Sevofluorane');
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

    fprintf('Data size: %dx%d (samples x channels)\n', size(data, 1), size(data, 2));
    fprintf('Fs = %.6f Hz, T = %.3f s\n', Fs, T);

    requestedElectrodeLabel = string(cfg.electrodeLabel);
    fallbackElectrodeLabel = "F4";

    f3Idx = find(strcmpi(labels, requestedElectrodeLabel), 1, 'first');
    f4Idx = find(strcmpi(labels, fallbackElectrodeLabel), 1, 'first');

    usedFallbackElectrode = false;
    electrodeSelectionNote = "";

    if ~isempty(f3Idx)
        if ~f3QualityPass
            error('F3 is present but its quality flag is false for %s. Exclude this subject from MFDB.', subjectID);
        end
        usedElectrodeLabel = requestedElectrodeLabel;
        usedElectrodeIndex = f3Idx;
        electrodeSelectionNote = "Used F3.";
    elseif ~isempty(f4Idx)
        usedElectrodeLabel = fallbackElectrodeLabel;
        usedElectrodeIndex = f4Idx;
        usedFallbackElectrode = true;
        electrodeSelectionNote = "F3 unavailable. Used F4 fallback.";
        warning('run_subject_mfdb_bootstrap:F4Fallback', ...
            'F3 not found for %s. Using F4 instead.', subjectID);
    else
        error('Neither F3 nor F4 was found for %s.', subjectID);
    end

    signalUsedElectrode = data(:, usedElectrodeIndex);

    expectedN = round(cfg.epochDurationSec * Fs);
    if N ~= expectedN
        warning('Expected N = %d for a %.1f s epoch, but found N = %d.', ...
            expectedN, cfg.epochDurationSec, N);
    end

    groupLabel = infer_group_label(subjectID);
    ageYears = lookup_subject_age_years(subjectID, cfg.anesthesia);

    fprintf('Requested electrode: %s\n', requestedElectrodeLabel);
    fprintf('Using %s index %d\n', usedElectrodeLabel, usedElectrodeIndex);
    if usedFallbackElectrode
        fprintf('Electrode note: %s\n', electrodeSelectionNote);
    end
    fprintf('Group: %s | AgeYears: %.4f\n', groupLabel, ageYears);

    [S_mtm, F_mtm, X_mtm, ~] = mtm(signalUsedElectrode, cfg.TW, cfg.K, Fs, cfg.method);

    nfft = size(X_mtm, 1);
    if mod(nfft, 2) ~= 0
        error('nfft must be even. Got %d.', nfft);
    end

    positiveIdx = (1:(nfft/2))';
    freqFull = F_mtm(positiveIdx);
    freqMask = freqFull >= cfg.fpass(1) & freqFull <= cfg.fpass(2);
    freq = freqFull(freqMask);
    nFreq = numel(freq);

    if nFreq == 0
        error('No frequency bins remain after fpass [%g %g].', cfg.fpass(1), cfg.fpass(2));
    end

    S_subject_original_linear = S_mtm(positiveIdx, 1);
    S_subject_original_linear = S_subject_original_linear(freqMask).';

    s = zeros(nfft, cfg.K);
    X_boot = complex(zeros(nfft/2, cfg.K, cfg.Bsubj));

    scale = sqrt(S_mtm(positiveIdx, 1) / 2);
    scale(scale == 0) = eps;

    s(1:(nfft/2), :) = real(X_mtm(positiveIdx, :, 1)) ./ scale;
    s((nfft/2 + 1):nfft, :) = imag(X_mtm(positiveIdx, :, 1)) ./ scale;

    s_mean = mean(s, 1);
    s_std = sqrt(mean((s - s_mean).^2, 1));
    s_std(s_std == 0) = eps;
    s_norm = (s - s_mean) ./ s_std;

    bootsam = randi(nfft, [nfft, cfg.Bsubj]);

    for k = 1:cfg.K
        s_boot = s_norm(:, k);
        s_boot = s_boot(bootsam);
        X_boot(:, k, :) = s_boot(1:(nfft/2), :) .* scale + ...
            1i * s_boot((nfft/2 + 1):nfft, :) .* scale;
    end

    PSD_boot = squeeze(mean(abs(X_boot).^2, 2));
    S_subject_boot_linear = PSD_boot(freqMask, :).';

    S_subject_original_dB = 10 * log10(max(S_subject_original_linear, eps));
    S_subject_boot_dB = 10 * log10(max(S_subject_boot_linear, eps));

    subjectOutDir = fullfile(cfg.outputDir, subjectID);
    if ~exist(subjectOutDir, 'dir')
        mkdir(subjectOutDir);
    end

    outFile = fullfile(subjectOutDir, sprintf('%s_%s_mfdb.mat', subjectID, upper(char(usedElectrodeLabel))));

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
    metadata.TW = cfg.TW;
    metadata.K = cfg.K;
    metadata.Bsubj = cfg.Bsubj;
    metadata.Bgroup = cfg.Bgroup;
    metadata.fpass = cfg.fpass;
    metadata.method = cfg.method;

    save(outFile, ...
        'subjectID', 'groupLabel', 'ageYears', 'freq', ...
        'S_subject_original_linear', 'S_subject_boot_linear', ...
        'S_subject_original_dB', 'S_subject_boot_dB', ...
        'Fs', 'T', 'N', 'requestedElectrodeLabel', 'usedElectrodeLabel', ...
        'usedElectrodeIndex', 'usedFallbackElectrode', 'electrodeSelectionNote', ...
        'metadata', '-v7.3');

    fprintf('Saved subject-level MFDB output to %s\n', outFile);

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
