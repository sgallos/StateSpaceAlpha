function plot_approach1_path_comparison()
%PLOT_APPROACH1_PATH_COMPARISON Four-panel comparison for Approach 1.
%
% Compares:
%   Path A, grid q      Path A, fixed q
%   Path B, grid q      Path B, fixed q
%
% All panels use MFDB observation noise only:
%   sigma2_obs = mfdb_var
%   sigma2_bio = 0

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outputsDir = fullfile(repoRoot, 'outputs');
figuresDir = fullfile(repoRoot, '4_figures');

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

variants = struct( ...
    'file', { ...
        'approach1_PathA_120s_results.mat', ...
        'approach1_PathA_120s_fixedq_results.mat', ...
        'approach1_PathB_subepoch_results.mat', ...
        'approach1_PathB_subepoch_fixedq_results.mat'}, ...
    'title', { ...
        'Path A, grid q', ...
        'Path A, fixed q', ...
        'Path B, grid q', ...
        'Path B, fixed q'}, ...
    'bandColor', { ...
        [0.60 0.72 0.92], ...
        [0.67 0.84 0.70], ...
        [0.93 0.74 0.58], ...
        [0.76 0.78 0.58]});

fig = figure('Color', 'w', 'Position', [50 50 1300 820]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for v = 1:numel(variants)
    resultFile = fullfile(outputsDir, variants(v).file);
    ax = nexttile(tl);
    hold(ax, 'on');

    if ~isfile(resultFile)
        text(ax, 0.5, 0.5, sprintf('Missing result file:\n%s', variants(v).file), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
        title(ax, variants(v).title, 'FontWeight', 'normal');
        axis(ax, 'off');
        continue;
    end

    S = load(resultFile, 'traj');
    traj = sortrows(S.traj, 'ageYears');
    ages = traj.ageYears;
    deltaMean = traj.deltaMean_dB;
    deltaLow = deltaMean - 1.96 * traj.deltaSD_dB;
    deltaHigh = deltaMean + 1.96 * traj.deltaSD_dB;

    patch(ax, [ages; flipud(ages)], [deltaLow; flipud(deltaHigh)], ...
        variants(v).bandColor, 'EdgeColor', 'none', 'FaceAlpha', 0.48);
    plot(ax, ages, deltaMean, '-k', 'LineWidth', 2);
    yline(ax, 0, '--', 'Color', [0.45 0.45 0.45]);

    title(ax, variants(v).title, 'FontWeight', 'normal');
    xlabel(ax, 'Age (years)');
    ylabel(ax, '\delta(a) (dB)');
    grid(ax, 'on');
    box(ax, 'off');
end

title(tl, 'Approach 1: MFDB observation noise only, grid-q vs fixed-q', ...
    'FontSize', 13, 'FontWeight', 'normal');

outFile = fullfile(figuresDir, 'approach1_four_panel_comparison.png');
exportgraphics(fig, outFile, 'Resolution', 200);
fprintf('Saved %s\n', outFile);
end
