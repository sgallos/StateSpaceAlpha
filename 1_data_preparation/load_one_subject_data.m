%% load_one_subject_data
% Minimal raw-data loader for one subject.
%
% After this runs, the workspace contains:
%   data      raw EEG data from the subject file
%   HDR       header struct from the subject file
%   dataPath  full path to the loaded file

clear; clc;

%% Pick one subject
subjectID = "CS29";

%% Find that subject's raw preprocessed file from the MFDB metadata
repoRoot = fileparts(fileparts(mfilename('fullpath')));
registryFile = fullfile(repoRoot, 'references', 'core_inputs', 'mfdb_subject_registry.csv');

registry = readtable(registryFile, 'TextType', 'string');
subjectRow = registry(strcmpi(registry.subjectID, subjectID), :);

if isempty(subjectRow)
    error('Subject %s was not found in %s.', subjectID, registryFile);
end

mfdbPath = subjectRow.filePath(1);
if ~isfile(mfdbPath)
    error('MFDB file not found for %s:\n%s', subjectID, mfdbPath);
end

mfdbInfo = load(mfdbPath, 'metadata');
if ~isfield(mfdbInfo, 'metadata') || ~isfield(mfdbInfo.metadata, 'dataPath')
    error('MFDB file for %s does not contain metadata.dataPath.', subjectID);
end

dataPath = mfdbInfo.metadata.dataPath;
if ~isfile(dataPath)
    error('Raw preprocessed file not found for %s:\n%s', subjectID, dataPath);
end

%% Load only the raw data and header
load(dataPath, 'data', 'HDR');

fprintf('Loaded %s\n', subjectID);
fprintf('File: %s\n', dataPath);
fprintf('data size: %s\n', mat2str(size(data)));
fprintf('HDR fields: %s\n', strjoin(string(fieldnames(HDR)).', ', '));
