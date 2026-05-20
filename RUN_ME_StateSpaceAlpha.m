function RUN_ME_StateSpaceAlpha()
%% RUN_ME_StateSpaceAlpha
% Main wrapper for the StateSpaceAlpha analysis.
%
% This is the human-facing entry point. It keeps the analysis run order in
% one place while leaving the actual analysis code in the numbered folders.
%
% Recommended meeting run:
%   1. Primary no-resampling posterior
%   2. Six-run no-resampling matrix
%   3. Refresh curated figures and CSVs

%% Section 1: Choose What To Run
% Flip these switches as needed. Defaults are the recommended current run.

runDataSmokeTest = false;              % Requires full subject MFDB files.
runCollapseSubepochTable = false;      % Use when the 4 x 30 s alpha table changed.

runPrimaryPosteriorNoResampling = true;
runNoResamplingMatrix = true;

% Diagnostics and secondary analyses. Leave off unless specifically needed.
runResidualSanityCheck = false;
runStep3LegacyNoEM = false;
runStep4JointEM = false;
runStep4FixedSigmaEM = false;
runStep4Subepoch76Diagnostic = false;
runStep5BootstrapSensitivity = false;

refreshCuratedReviewFiles = true;
makeFiguresVisible = true;            % true shows figures while scripts run.

%% Section 2: Set Up Paths
repoRoot = fileparts(mfilename('fullpath'));
outputsDir = fullfile(repoRoot, 'outputs');
figuresDir = fullfile(repoRoot, '4_figures');
curatedDataDir = fullfile(repoRoot, '5_outputs_data');

addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));

if ~exist(outputsDir, 'dir'), mkdir(outputsDir); end
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end
if ~exist(curatedDataDir, 'dir'), mkdir(curatedDataDir); end

if makeFiguresVisible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

fprintf('\nStateSpaceAlpha main wrapper\n');
fprintf('Repository: %s\n', repoRoot);
fprintf('Figures visible: %d\n', makeFiguresVisible);

%% Section 3: Data Preparation
if runDataSmokeTest
    run_child_script(repoRoot, '1_data_preparation', 'smoke_test_alpha_extraction.m');
else
    fprintf('\nSkipping smoke-test alpha extraction.\n');
end

if runCollapseSubepochTable
    run_child_script(repoRoot, '1_data_preparation', 'collapse_subepoch_alpha_table_for_ssm.m');
else
    fprintf('Skipping sub-epoch table collapse.\n');
end

%% Section 4: Primary No-Resampling SSM Analyses
if runPrimaryPosteriorNoResampling
    run_child_script(repoRoot, '3_run_analyses', 'run_ssm_posterior_no_resampling.m');
else
    fprintf('\nSkipping primary no-resampling posterior workflow.\n');
end

if runNoResamplingMatrix
    run_child_script(repoRoot, '3_run_analyses', 'run_ssm_no_resampling_matrix.m');
else
    fprintf('\nSkipping six-run no-resampling matrix.\n');
end

%% Section 5: Optional Diagnostics
if runResidualSanityCheck
    run_child_script(repoRoot, '3_run_analyses', 'run_ssm_residual_sanity_check.m');
end

if runStep3LegacyNoEM
    run_child_script(repoRoot, '3_run_analyses', 'run_ss_age_diff.m');
end

if runStep4JointEM
    run_child_script(repoRoot, '3_run_analyses', 'run_ss_age_diff_em.m');
end

if runStep4FixedSigmaEM
    run_child_script(repoRoot, '3_run_analyses', 'run_ss_age_diff_em_fixed_sigma.m');
end

if runStep4Subepoch76Diagnostic
    run_child_script(repoRoot, '3_run_analyses', 'run_ss_age_diff_em_subepoch76.m');
end

if runStep5BootstrapSensitivity
    run_child_script(repoRoot, '3_run_analyses', 'run_ss_age_diff_bootstrap_simultaneous.m');
end

%% Section 6: Refresh Curated Review Files
if refreshCuratedReviewFiles
    refresh_curated_review_files(outputsDir, figuresDir, curatedDataDir);
else
    fprintf('\nSkipping curated review-file refresh.\n');
end

%% Section 7: Final Pointers
fprintf('\nStateSpaceAlpha wrapper complete.\n');
fprintf('Read this first: %s\n', fullfile(repoRoot, 'RESULTS_SUMMARY.md'));
fprintf('Curated figures: %s\n', figuresDir);
fprintf('Curated CSVs: %s\n', curatedDataDir);

end

function run_child_script(repoRoot, folderName, scriptName)
    scriptPath = fullfile(repoRoot, folderName, scriptName);
    if ~isfile(scriptPath)
        error('Cannot find script: %s', scriptPath);
    end

    fprintf('\n============================================================\n');
    fprintf('Running %s/%s\n', folderName, scriptName);
    fprintf('============================================================\n');

    % Child scripts start with clear; clc. Running them in base keeps those
    % commands from clearing this wrapper function's settings.
    evalin('base', sprintf('run(''%s'')', escape_for_matlab_string(scriptPath)));
end

function refresh_curated_review_files(outputsDir, figuresDir, curatedDataDir)
    fprintf('\nRefreshing curated review files...\n');

    copy_if_present(outputsDir, figuresDir, ...
        'ssm_posterior_delta_primary_no_resampling.png', ...
        'primary_delta_trajectory.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_no_resampling_matrix_delta_overlay.png', ...
        'no_resampling_matrix_overlay.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_no_resampling_matrix_delta_bands.png', ...
        'no_resampling_matrix_delta_bands.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_no_resampling_matrix_profile_loglik_heatmap.png', ...
        'profile_likelihood_heatmap.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_step4_em_loglik_history.png', ...
        'em_loglik_history.png');

    copy_if_present(outputsDir, curatedDataDir, ...
        'ssm_no_resampling_matrix_summary.csv', ...
        'ssm_no_resampling_matrix_summary.csv');
    copy_if_present(outputsDir, curatedDataDir, ...
        'ssm_no_resampling_matrix_trajectories.csv', ...
        'ssm_no_resampling_matrix_trajectories.csv');
    copy_if_present(outputsDir, curatedDataDir, ...
        'ssm_posterior_no_resampling_primary_trajectory.csv', ...
        'ssm_posterior_no_resampling_primary_trajectory.csv');
    copy_if_present(outputsDir, curatedDataDir, ...
        'ssm_posterior_no_resampling_sensitivity_summary.csv', ...
        'ssm_posterior_no_resampling_sensitivity_summary.csv');
    copy_if_present(outputsDir, curatedDataDir, ...
        'ssm_posterior_no_resampling_sigmaBio_subject_diagnostics.csv', ...
        'ssm_posterior_no_resampling_sigmaBio_subject_diagnostics.csv');
end

function copy_if_present(sourceDir, destinationDir, sourceName, destinationName)
    sourcePath = fullfile(sourceDir, sourceName);
    destinationPath = fullfile(destinationDir, destinationName);

    if isfile(sourcePath)
        copyfile(sourcePath, destinationPath);
        fprintf('  copied %s -> %s\n', sourceName, destinationName);
    else
        fprintf('  missing %s; skipped\n', sourceName);
    end
end

function escaped = escape_for_matlab_string(pathText)
    escaped = strrep(pathText, '''', '''''');
end
