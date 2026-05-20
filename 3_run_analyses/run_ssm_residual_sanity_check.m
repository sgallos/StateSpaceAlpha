%% run_ssm_residual_sanity_check
% Quick diagnostic for whether MFDB observation variance is much smaller
% than residual scatter around the current Step 3 fitted trajectories.

clear; clc;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(repoRoot, 'outputs');
addpath(fullfile(repoRoot, '1_data_preparation'));
addpath(fullfile(repoRoot, '2_state_space_model'));
addpath(fullfile(repoRoot, '3_run_analyses'));
step3MatFile = fullfile(outDir, 'ssm_step3_age_difference_results.mat');

if ~isfile(step3MatFile)
    error('Step 3 results not found. Run run_ss_age_diff first.');
end

S = load(step3MatFile, 'mainResult');
data = S.mainResult.subjectTable;
trajectory = S.mainResult.trajectory;

fittedAlpha = trajectory.baselineMean_dB + ...
    trajectory.groupIndicator .* trajectory.deltaMean_dB;
residuals = data.alpha_dB - fittedAlpha;

residualSD = std(residuals);
meanSqrtMfdb = mean(sqrt(data.mfdb_var));
impliedBiologicalVariance = max(residualSD.^2 - meanSqrtMfdb.^2, 0);
impliedBiologicalSD = sqrt(impliedBiologicalVariance);

fprintf('\nResidual SD from current fit: %.6f dB\n', residualSD);
fprintf('Mean sqrt(mfdb_var):          %.6f dB\n', meanSqrtMfdb);
fprintf('Implied biological variance:  %.6f dB^2\n', impliedBiologicalVariance);
fprintf('Implied biological SD:        %.6f dB\n', impliedBiologicalSD);

fig = figure('Color', 'w', 'Position', [100 100 1100 520]);
tiledlayout(fig, 1, 2);

ax1 = nexttile;
histogram(ax1, data.alpha_dB, 10);
xlabel(ax1, 'Alpha power (dB)');
ylabel(ax1, 'Count');
title(ax1, 'Alpha power distribution');

ax2 = nexttile;
gscatter(data.ageYears, data.alpha_dB, data.groupLabel);
hold(ax2, 'on');
plot(ax2, trajectory.ageYears, trajectory.baselineMean_dB, ...
    '-k', 'LineWidth', 2);
plot(ax2, trajectory.ageYears, ...
    trajectory.baselineMean_dB + trajectory.deltaMean_dB, ...
    '--k', 'LineWidth', 2);
xlabel(ax2, 'Age (years)');
ylabel(ax2, 'Alpha power (dB)');
title(ax2, 'Alpha power vs age with Step 3 fitted curves');
legend(ax2, 'Location', 'best');

sanityFigureFile = fullfile(outDir, 'ssm_residual_sanity_check.png');
saveas(fig, sanityFigureFile);
fprintf('Saved sanity-check figure to %s\n', sanityFigureFile);
