function plot_em_convergence_diagnostics(matrixResultsFile, outputFile)
%PLOT_EM_CONVERGENCE_DIAGNOSTICS Visualize EM behavior across strategies.
%
% This diagnostic reads the no-resampling matrix MAT file and plots the four
% EM-based strategies:
%   EM-Smoothness
%   EM-Joint
%   EM-Joint-Subepoch
%   EM-Biological
%
% Rows:
%   1. log-likelihood history
%   2. q_f0 and q_delta histories
%   3. sigmaBio history

if nargin < 1 || isempty(matrixResultsFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    matrixResultsFile = fullfile(repoRoot, 'outputs', 'ssm_no_resampling_matrix_results.mat');
end

if nargin < 2 || isempty(outputFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    outputFile = fullfile(repoRoot, 'outputs', 'em_convergence_diagnostics.png');
end

if ~isfile(matrixResultsFile)
    error(['Missing no-resampling matrix results file:\n%s\n\n' ...
        'Run 3_run_analyses/run_ssm_no_resampling_matrix.m first.'], matrixResultsFile);
end

S = load(matrixResultsFile);

strategyList = struct( ...
    'displayName', {'EM-Smoothness', 'EM-Joint', 'EM-Joint-Subepoch', 'EM-Biological'}, ...
    'resultField', {'emSmoothnessResults', 'emJointResults', 'emJointSubepochResults', 'emBiologicalResults'});
nStrategies = numel(strategyList);

fig = figure('Color', 'w', 'Position', [50 50 1550 870]);
layout = tiledlayout(fig, 3, nStrategies, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:nStrategies
    strategy = strategyList(s);
    emResult = extract_em_result(S, strategy.resultField);

    if isempty(emResult)
        plot_missing_strategy_column(strategy.displayName, s, nStrategies);
        continue;
    end

    history = standardize_em_history(emResult);
    status = summarize_em_status(emResult, history);

    plot_log_likelihood_panel(s, nStrategies, strategy.displayName, history, status);
    plot_q_history_panel(s, nStrategies, history, status);
    plot_sigma_bio_panel(s, nStrategies, history);
end

title(layout, 'EM convergence diagnostics across no-resampling strategies', ...
    'FontSize', 14, 'FontWeight', 'bold');

outputDir = fileparts(outputFile);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
exportgraphics(fig, outputFile, 'Resolution', 200);
fprintf('Saved EM convergence diagnostics to %s\n', outputFile);
end

function emResult = extract_em_result(S, resultField)
    emResult = [];
    if isfield(S, resultField)
        emResult = S.(resultField);
    end
end

function history = standardize_em_history(emResult)
    bestStart = emResult.bestStart;

    history.logLik = bestStart.logLikHistory(:);
    history.q = bestStart.qHistory;
    history.sigmaBio = bestStart.sigmaBioHistory(:);

    history.logLik = history.logLik(isfinite(history.logLik));
    history.q = history.q(all(isfinite(history.q), 2), :);
    history.sigmaBio = history.sigmaBio(isfinite(history.sigmaBio));

    nIter = max([numel(history.logLik), size(history.q, 1), numel(history.sigmaBio)]);
    history.iterations = (1:nIter).';
end

function status = summarize_em_status(emResult, history)
    status = struct();
    status.qCeiling = get_option_or_default(emResult, 'qCeiling', 1e6);
    status.logLikIsMonotone = all(diff(history.logLik) >= -1e-7);
    status.hitQCeiling = any(history.q(:) >= 0.999 * status.qCeiling) || ...
        safe_get_logical(emResult, 'hitQCeiling');
    status.hitQFloor = safe_get_logical(emResult, 'hitQFloor');
    status.hitSigmaBioCeiling = safe_get_logical(emResult, 'hitSigmaBioCeiling');
    status.hitSigmaBioFloor = safe_get_logical(emResult, 'hitSigmaBioFloor');
    status.hitBoundary = status.hitQCeiling || status.hitQFloor || ...
        status.hitSigmaBioCeiling || status.hitSigmaBioFloor;
    status.converged = safe_get_logical(emResult, 'converged');
    status.healthy = status.logLikIsMonotone && status.converged && ~status.hitBoundary;
end

function plot_log_likelihood_panel(tileIndex, nStrategies, displayName, history, status)
    ax = nexttile(tileIndex);
    plot(ax, 1:numel(history.logLik), history.logLik, '-o', ...
        'LineWidth', 1.4, 'MarkerSize', 4);

    if status.healthy
        statusText = 'OK';
        titleColor = [0 0.45 0];
    else
        statusText = 'ISSUE';
        titleColor = [0.72 0 0];
    end

    title(ax, sprintf('%s (%s)', displayName, statusText), ...
        'Color', titleColor, 'Interpreter', 'none');
    ylabel(ax, 'log-likelihood');

    statusLines = make_status_lines(status);
    text(ax, 0.03, 0.95, cellstr(statusLines), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 8, ...
        'Color', titleColor, 'Interpreter', 'none');

    style_axis(ax);
    if tileIndex == 1
        add_row_label(ax, 'Log-likelihood');
    end
    maybe_hide_x_tick_labels(ax, tileIndex, nStrategies);
end

function plot_q_history_panel(tileIndex, nStrategies, history, status)
    ax = nexttile(tileIndex + nStrategies);
    semilogy(ax, 1:size(history.q, 1), history.q(:, 1), '-o', ...
        'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'q_f0');
    hold(ax, 'on');
    semilogy(ax, 1:size(history.q, 1), history.q(:, 2), '-s', ...
        'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', 'q_delta');
    yline(ax, status.qCeiling, '--', 'ceiling', ...
        'Color', [0.72 0 0], 'HandleVisibility', 'off');
    ylabel(ax, 'q (log scale)');
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
    style_axis(ax);
    if tileIndex == 1
        add_row_label(ax, 'Smoothness');
    end
    maybe_hide_x_tick_labels(ax, tileIndex, nStrategies);
end

function plot_sigma_bio_panel(tileIndex, nStrategies, history)
    ax = nexttile(tileIndex + 2 * nStrategies);
    plot(ax, 1:numel(history.sigmaBio), history.sigmaBio, '-o', ...
        'LineWidth', 1.4, 'MarkerSize', 4);
    ylabel(ax, 'sigmaBio (dB^2)');
    xlabel(ax, 'EM iteration');

    if numel(history.sigmaBio) > 1 && all(abs(diff(history.sigmaBio)) < 1e-12)
        text(ax, 0.03, 0.95, 'fixed', 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 8, ...
            'Color', [0.3 0.3 0.3]);
    end

    style_axis(ax);
    if tileIndex == 1
        add_row_label(ax, 'Biology variance');
    end
end

function plot_missing_strategy_column(displayName, tileIndex, nStrategies)
    for row = 1:3
        ax = nexttile(tileIndex + (row - 1) * nStrategies);
        text(ax, 0.5, 0.5, sprintf('No data for\n%s', displayName), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'none');
        axis(ax, 'off');
    end
end

function lines = make_status_lines(status)
    lines = strings(0);
    if status.logLikIsMonotone
        lines(end + 1) = "log-lik monotonic";
    else
        lines(end + 1) = "log-lik non-monotonic";
    end

    if status.converged
        lines(end + 1) = "converged";
    else
        lines(end + 1) = "did not converge";
    end

    if status.hitQCeiling
        lines(end + 1) = "q hit ceiling";
    elseif status.hitQFloor
        lines(end + 1) = "q hit floor";
    end

    if status.hitSigmaBioCeiling
        lines(end + 1) = "sigmaBio hit ceiling";
    elseif status.hitSigmaBioFloor
        lines(end + 1) = "sigmaBio hit floor";
    end
end

function value = get_option_or_default(emResult, optionName, defaultValue)
    value = defaultValue;
    if isfield(emResult, 'options') && isfield(emResult.options, optionName)
        value = emResult.options.(optionName);
    end
end

function value = safe_get_logical(emResult, fieldName)
    value = false;
    if isfield(emResult, fieldName)
        value = logical(emResult.(fieldName));
    elseif isfield(emResult, 'bestStart') && isfield(emResult.bestStart, fieldName)
        value = logical(emResult.bestStart.(fieldName));
    end
end

function add_row_label(ax, labelText)
    text(ax, -0.22, 0.5, labelText, 'Units', 'normalized', ...
        'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
        'Interpreter', 'none');
end

function maybe_hide_x_tick_labels(ax, tileIndex, nStrategies)
    if tileIndex <= 2 * nStrategies
        ax.XTickLabel = [];
    end
end

function style_axis(ax)
    ax.FontName = 'Helvetica';
    ax.FontSize = 10;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    box(ax, 'off');
    grid(ax, 'on');
end
