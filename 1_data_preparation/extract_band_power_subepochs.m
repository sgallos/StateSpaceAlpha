function T = extract_band_power_subepochs()
%EXTRACT_BAND_POWER_SUBEPOCHS Build a 76-row delta/theta/alpha power table.
%
% Uses the saved 4 x 30 s sub-epoch MFDB spectra. No multitaper or MFDB
% bootstrap is rerun here; the script only reintegrates the saved spectra
% over each band:
%   delta: 1-4 Hz
%   theta: 4-8 Hz
%   alpha: 8-12 Hz
%
% Output:
%   outputs/power_table_for_ssm_subepochs.csv
%   outputs/power_table_for_ssm_subepochs.mat

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

multitaperRoot = '/Users/gallo/Documents/MATLAB_Multitaper_Hz_Domain_BTS';
oldPath = path;
cleanupObj = onCleanup(@() path(oldPath));
addpath(multitaperRoot);
cfg = mfdb_config();
subepochRoot = cfg.subEpochOutputDir;
validationFile = cfg.validationFile;

if ~isfile(validationFile)
    error('Validation file not found: %s', validationFile);
end
if ~isfolder(subepochRoot)
    error('Sub-epoch output folder not found: %s', subepochRoot);
end

bands = struct( ...
    'name', {'delta', 'theta', 'alpha'}, ...
    'lowHz', {1, 4, 8}, ...
    'highHz', {4, 8, 12});

V = load(validationFile, 'includedAfterValidation');
registry = V.includedAfterValidation;
nSubjects = height(registry);
nExpectedRows = nSubjects * 4;

T = table('Size', [nExpectedRows, 13], ...
    'VariableTypes', {'string','double','string','double', ...
    'double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'subjectID','subepochIdx','groupLabel','ageYears', ...
    'delta_dB','theta_dB','alpha_dB', ...
    'delta_mfdb_var','theta_mfdb_var','alpha_mfdb_var', ...
    'delta_nFreqBins','theta_nFreqBins','alpha_nFreqBins'});

rowIdx = 0;
for subjectIdx = 1:nSubjects
    subjectID = string(registry.subjectID(subjectIdx));
    subepochFile = find_subepoch_file(subepochRoot, subjectID);

    if subepochFile == ""
        warning('Missing sub-epoch file for %s. Skipping.', subjectID);
        continue;
    end

    S = load(subepochFile, 'subEpochResults', 'groupLabel', 'ageYears');
    for subepochIdx = 1:numel(S.subEpochResults)
        rowIdx = rowIdx + 1;
        thisSubEpoch = S.subEpochResults{subepochIdx};

        T.subjectID(rowIdx) = subjectID;
        T.subepochIdx(rowIdx) = subepochIdx;
        T.groupLabel(rowIdx) = string(S.groupLabel);
        T.ageYears(rowIdx) = double(S.ageYears);

        for bandIdx = 1:numel(bands)
            band = bands(bandIdx);
            [powerDB, mfdbVar, nFreqBins] = integrate_subepoch_band(thisSubEpoch, ...
                band.lowHz, band.highHz, subjectID, subepochIdx, band.name);

            bandName = string(band.name);
            T.(bandName + "_dB")(rowIdx) = powerDB;
            T.(bandName + "_mfdb_var")(rowIdx) = mfdbVar;
            T.(bandName + "_nFreqBins")(rowIdx) = nFreqBins;
        end
    end
end

T = T(1:rowIdx, :);

outCsv = fullfile(outDir, 'power_table_for_ssm_subepochs.csv');
outMat = fullfile(outDir, 'power_table_for_ssm_subepochs.mat');
writetable(T, outCsv);
save(outMat, 'T', 'bands', 'subepochRoot', 'validationFile');

fprintf('Wrote %d sub-epoch rows to %s\n', height(T), outCsv);
disp(T(1:min(12, height(T)), :));
end

function subepochFile = find_subepoch_file(subepochRoot, subjectID)
subjectDir = fullfile(subepochRoot, char(subjectID));
files = dir(fullfile(subjectDir, sprintf('%s_*_mfdb_subepochs.mat', char(subjectID))));

if isempty(files)
    subepochFile = "";
else
    subepochFile = string(fullfile(files(1).folder, files(1).name));
end
end

function [powerDB, mfdbVar, nFreqBins] = integrate_subepoch_band(thisSubEpoch, ...
    lowHz, highHz, subjectID, subepochIdx, bandName)
freq = double(thisSubEpoch.freq(:).');
bandMask = freq >= lowHz & freq <= highHz;
nFreqBins = sum(bandMask);

if nFreqBins == 0
    error('No %s bins found for %s sub-epoch %d in [%g, %g] Hz.', ...
        bandName, subjectID, subepochIdx, lowHz, highHz);
end

originalDB = double(thisSubEpoch.S_subject_original_dB);
bootDB = double(thisSubEpoch.S_subject_boot_dB);

originalLinear = trapz(freq(bandMask), 10.^(originalDB(bandMask) / 10));
powerDB = 10 * log10(max(originalLinear, eps));

bootLinear = trapz(freq(bandMask), 10.^(bootDB(:, bandMask) / 10), 2);
bootPowerDB = 10 * log10(max(bootLinear, eps));
mfdbVar = var(bootPowerDB, 0, 'omitnan');
end
