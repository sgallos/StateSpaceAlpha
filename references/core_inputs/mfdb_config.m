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
end
