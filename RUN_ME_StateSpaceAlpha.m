function RUN_ME_StateSpaceAlpha()
%% RUN_ME_StateSpaceAlpha
% Full analysis pipeline for the CP-vs-Control alpha power study.
%
% The analysis tells four chapters:
%
%   Chapter 1: Get the data ready.
%       Goal: turn MFDB outputs into one clean per-subject alpha table.
%
%   Chapter 2: Look at the data before modeling.
%       Goal: see what the model sees before imposing a smoothness prior.
%
%   Chapter 3: Fit the primary state-space model.
%       Goal: produce the CP-minus-Control alpha trajectory with credible bands.
%
%   Chapter 4: Check sensitivity to model choices.
%       Goal: show whether the finding is stable across six hyperparameter strategies.

%% ====================================================================
%% CHOOSE WHICH CHAPTERS TO RUN
%% ====================================================================
% Defaults are set for a normal collaborator-facing rerun.

runChapter1_prepareData = false;        % Rebuilds alpha tables; slow if upstream files are missing.
runChapter2_plotRawData = true;         % Subject-level alpha-vs-age scatter.
runChapter3_primaryModel = true;        % Primary SSM posterior fit.
runChapter4_sensitivityMatrix = true;   % Six-strategy robustness matrix.

% Optional diagnostics and secondary analyses. Leave off unless needed.
runResidualSanityCheck = false;
runStep3LegacyNoEM = false;
runStep4JointEM = false;
runStep4FixedSigmaEM = false;
runStep4Subepoch76Diagnostic = false;
runStep5BootstrapSensitivity = false;

refreshCuratedReviewFiles = true;
makeFiguresVisible = true;              % true shows figures while scripts run.

%% ====================================================================
%% SET UP PATHS
%% ====================================================================
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

fprintf('\nStateSpaceAlpha full analysis pipeline\n');
fprintf('Repository: %s\n', repoRoot);
fprintf('Figures visible: %d\n', makeFiguresVisible);

%% ====================================================================
%% CHAPTER 1: GET THE DATA READY
%% ====================================================================
fprintf('\n=== CHAPTER 1: Prepare the alpha table ===\n');

if runChapter1_prepareData
    % 1.1 Extract alpha power from the existing 120 s subject-level MFDB
    % outputs. This is a smoke test that confirms our alpha integration
    % matches the original pipeline.
    %
    % Produces:
    %   outputs/alpha_table_for_ssm.csv
    %   outputs/smoke_test_alpha_vs_age.png
    run_child_script(repoRoot, '1_data_preparation', 'smoke_test_alpha_extraction.m');

    % 1.2 Check whether each subject's four 30 s alpha estimates are
    % internally consistent.
    %
    % Produces:
    %   outputs/subepoch_consistency_summary.csv
    %   outputs/subepoch_consistency_scatter.png
    run_child_script(repoRoot, '1_data_preparation', 'check_subepoch_consistency.m');

    % 1.3 Collapse the 4 x 30 s sub-epoch alpha table into one row per
    % subject. The collapsed mfdb_var uses sub-epoch variability, which is
    % more realistic than the original 120 s MFDB-only variance.
    %
    % Produces:
    %   outputs/alpha_table_for_ssm_subepochs_collapsed.csv
    %   outputs/alpha_table_subepoch_vs_120s_comparison.png
    run_child_script(repoRoot, '1_data_preparation', 'collapse_subepoch_alpha_table_for_ssm.m');
else
    fprintf('Skipping Chapter 1. Existing alpha tables will be used.\n');
end

%% ====================================================================
%% CHAPTER 2: LOOK AT THE DATA BEFORE MODELING
%% ====================================================================
fprintf('\n=== CHAPTER 2: Eyeball the data ===\n');

if runChapter2_plotRawData
    % 2.1 Plot subject-level alpha vs age, colored by group.
    %
    % Figure:
    %   outputs/subject_level_alpha_vs_age.png
    %
    % What to look for:
    %   19 points, with controls and CP separated by color. This is the
    %   actual collapsed subject-level table the SSM sees.
    run_child_script(repoRoot, '3_run_analyses', 'plot_subject_level_alpha_vs_age.m');
else
    fprintf('Skipping Chapter 2 raw-data plot.\n');
end

%% ====================================================================
%% CHAPTER 3: PRIMARY STATE-SPACE MODEL FIT
%% ====================================================================
fprintf('\n=== CHAPTER 3: Primary SSM fit ===\n');

if runChapter3_primaryModel
    % 3.1 Fit the SSM with empirical biological variance and fixed
    % smoothness. No EM, no subject resampling. Credible bands come directly
    % from the smoother posterior covariance.
    %
    % Figure:
    %   outputs/ssm_posterior_delta_primary_no_resampling.png
    %
    % Key interpretation:
    %   The headline CP-minus-Control trajectory. Non-overlap with zero is
    %   where CP differs detectably from Control under this model.
    run_child_script(repoRoot, '3_run_analyses', 'run_ssm_posterior_no_resampling.m');

    % 3.2 Plot the actual Control and CP posterior trajectories with the
    % 19 subject-level alpha estimates overlaid. This is the concrete data
    % view that complements the formal delta(a) plot.
    %
    % Figure:
    %   outputs/subjects_with_trajectories.png
    plot_subjects_with_trajectories( ...
        fullfile(outputsDir, 'ssm_posterior_no_resampling_results.mat'), ...
        fullfile(outputsDir, 'subjects_with_trajectories.png'));
else
    fprintf('Skipping Chapter 3 primary SSM fit.\n');
end

%% ====================================================================
%% CHAPTER 4: SENSITIVITY TO HYPERPARAMETER CHOICE
%% ====================================================================
fprintf('\n=== CHAPTER 4: Sensitivity matrix ===\n');

if runChapter4_sensitivityMatrix
    % 4.1 Run six hyperparameter selection strategies. All use the SSM
    % smoother posterior. No subject resampling.
    %
    % Strategies:
    %   Fixed-Heuristic     - all hyperparameters at sensible defaults
    %   EM-Smoothness       - EM finds smoothness, biology fixed empirically
    %   EM-Joint            - EM finds everything jointly
    %   EM-Joint-Subepoch   - same idea on the 76-row table
    %   EM-Biological       - EM finds biology, smoothness fixed
    %   Profile-Likelihood  - grid search for maximum-likelihood smoothness
    %
    % Figures:
    %   outputs/sensitivity_all_six_overlay.png
    %   outputs/sensitivity_all_six_bands.png
    %   outputs/sensitivity_profile_loglik_heatmap.png
    %   outputs/em_convergence_diagnostics.png
    run_child_script(repoRoot, '3_run_analyses', 'run_ssm_no_resampling_matrix.m');

    % 4.2 Plot the EM histories behind the EM-based sensitivity strategies.
    % This shows why some strategies are clean and others are diagnostic
    % failures at n = 19.
    plot_em_convergence_diagnostics( ...
        fullfile(outputsDir, 'ssm_no_resampling_matrix_results.mat'), ...
        fullfile(outputsDir, 'em_convergence_diagnostics.png'));
else
    fprintf('Skipping Chapter 4 sensitivity matrix.\n');
end

%% ====================================================================
%% OPTIONAL DIAGNOSTICS
%% ====================================================================
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

%% ====================================================================
%% REFRESH CURATED REVIEW FILES
%% ====================================================================
if refreshCuratedReviewFiles
    refresh_curated_review_files(outputsDir, figuresDir, curatedDataDir);
else
    fprintf('\nSkipping curated review-file refresh.\n');
end

%% ====================================================================
%% DONE
%% ====================================================================
fprintf('\nDone.\n');
fprintf('Read RESULTS_SUMMARY.md for the writeup.\n');
fprintf('Browse 4_figures/ for the plots.\n');
fprintf('See 4_figures/README.md for what each figure shows.\n');

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
        'subject_level_alpha_vs_age.png', ...
        'subject_level_alpha_vs_age.png');
    copy_if_present(outputsDir, figuresDir, ...
        'alpha_table_subepoch_vs_120s_comparison.png', ...
        'alpha_table_subepoch_vs_120s_comparison.png');
    copy_if_present(outputsDir, figuresDir, ...
        'subepoch_consistency_scatter.png', ...
        'subepoch_consistency_scatter.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_posterior_delta_primary_no_resampling.png', ...
        'primary_delta_trajectory.png');
    copy_if_present(outputsDir, figuresDir, ...
        'subjects_with_trajectories.png', ...
        'subjects_with_trajectories.png');
    copy_if_present(outputsDir, figuresDir, ...
        'sensitivity_all_six_overlay.png', ...
        'sensitivity_all_six_overlay.png');
    copy_if_present(outputsDir, figuresDir, ...
        'sensitivity_all_six_bands.png', ...
        'sensitivity_all_six_bands.png');
    copy_if_present(outputsDir, figuresDir, ...
        'sensitivity_profile_loglik_heatmap.png', ...
        'sensitivity_profile_loglik_heatmap.png');
    copy_if_present(outputsDir, figuresDir, ...
        'em_convergence_diagnostics.png', ...
        'em_convergence_diagnostics.png');
    copy_if_present(outputsDir, figuresDir, ...
        'ssm_step4_em_loglik_history.png', ...
        'em_loglik_history.png');

    copy_if_present(outputsDir, curatedDataDir, ...
        'subepoch_consistency_summary.csv', ...
        'subepoch_consistency_summary.csv');
    copy_if_present(outputsDir, curatedDataDir, ...
        'alpha_table_subepoch_vs_120s_comparison.csv', ...
        'alpha_table_subepoch_vs_120s_comparison.csv');
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
