function cfg = mfdb_config()
%MFDB_CONFIG Locked configuration for subject- and group-level MFDB analysis.

    cfg = struct();

    cfg.scriptRoot = '/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS';
    cfg.pedsRoot = '/Users/gallo/Documents/peds_cp';

    cfg.anesthesia = 'sevo';
    cfg.electrodeLabel = 'F3';
    cfg.epochDurationSec = 120;
    cfg.TW = 6;
    cfg.K = 11;
    cfg.Bsubj = 2000;
    cfg.Bgroup = 2000;
    cfg.fpass = [1 40];
    cfg.method = 'unity';

    cfg.outputDir = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_subject_level');
    cfg.validationFile = fullfile(cfg.scriptRoot, 'outputs', 'mfdb_validation', ...
        'mfdb_validation_results.mat');

    % Sub-epoch MFDB settings. These are used only by the parallel
    % sub-epoch pipeline and do not change the original 120 s subject files.
    cfg.subEpochDurationSec = 30;
    cfg.subEpochCount = 4;
    cfg.subEpochTW = 4;
    cfg.subEpochK = 7;
    cfg.subEpochOutputDir = fullfile(cfg.scriptRoot, 'outputs', ...
        'mfdb_subepoch_level');

    assert(cfg.subEpochCount * cfg.subEpochDurationSec == cfg.epochDurationSec, ...
        'Sub-epoch count and duration must sum to the original epoch duration.');
end
