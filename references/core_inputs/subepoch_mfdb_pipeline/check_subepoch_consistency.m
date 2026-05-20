function perSubject = check_subepoch_consistency()
%CHECK_SUBEPOCH_CONSISTENCY Summarize within-subject alpha variability.
%
% This is the main QC step for deciding whether four 30 s sub-epochs are
% stable enough to use as repeated observations.

    cfg = mfdb_config();
    alphaTableFile = fullfile(cfg.scriptRoot, 'outputs', ...
        'alpha_table_for_ssm_subepochs.csv');

    if ~isfile(alphaTableFile)
        error(['Sub-epoch alpha table not found:\n%s\n\n' ...
            'Run smoke_test_alpha_extraction_subepochs first.'], alphaTableFile);
    end

    T = readtable(alphaTableFile, 'TextType', 'string');
    subjectIDs = unique(T.subjectID, 'stable');
    nSubjects = numel(subjectIDs);

    perSubject = table();
    perSubject.subjectID = strings(nSubjects, 1);
    perSubject.groupLabel = strings(nSubjects, 1);
    perSubject.ageYears = nan(nSubjects, 1);
    perSubject.nSubEpochs = zeros(nSubjects, 1);
    perSubject.alpha_mean = nan(nSubjects, 1);
    perSubject.alpha_sd = nan(nSubjects, 1);
    perSubject.alpha_range = nan(nSubjects, 1);
    perSubject.mean_mfdb_sd = nan(nSubjects, 1);

    for i = 1:nSubjects
        rows = T.subjectID == subjectIDs(i);
        thisData = sortrows(T(rows, :), 'subepochIdx');
        alphaVals = thisData.alpha_dB;

        perSubject.subjectID(i) = subjectIDs(i);
        perSubject.groupLabel(i) = thisData.groupLabel(1);
        perSubject.ageYears(i) = thisData.ageYears(1);
        perSubject.nSubEpochs(i) = height(thisData);
        perSubject.alpha_mean(i) = mean(alphaVals, 'omitnan');
        perSubject.alpha_sd(i) = std(alphaVals, 'omitnan');
        perSubject.alpha_range(i) = max(alphaVals) - min(alphaVals);
        perSubject.mean_mfdb_sd(i) = mean(sqrt(thisData.mfdb_var), 'omitnan');
    end

    perSubject = sortrows(perSubject, 'alpha_sd', 'descend');

    outDir = fullfile(cfg.scriptRoot, 'outputs');
    summaryFile = fullfile(outDir, 'subepoch_consistency_summary.csv');
    plotFile = fullfile(outDir, 'subepoch_consistency_scatter.png');

    writetable(perSubject, summaryFile);

    fprintf('Per-subject sub-epoch consistency:\n');
    disp(perSubject);
    fprintf('\nMedian within-subject SD: %.3f dB\n', ...
        median(perSubject.alpha_sd, 'omitnan'));
    fprintf('Median within-subject range: %.3f dB\n', ...
        median(perSubject.alpha_range, 'omitnan'));
    fprintf('Median mean MFDB SD: %.4f dB\n', ...
        median(perSubject.mean_mfdb_sd, 'omitnan'));
    fprintf('Saved consistency summary to %s\n', summaryFile);

    make_consistency_plot(T, subjectIDs, plotFile);
    fprintf('Saved consistency scatter plot to %s\n', plotFile);
end

function make_consistency_plot(T, subjectIDs, plotFile)
    fig = figure('Color', 'w', 'Position', [100 100 1050 620]);
    ax = axes(fig);
    hold(ax, 'on');

    controlColor = [0.00 0.45 0.74];
    cpColor = [0.85 0.33 0.10];

    for i = 1:numel(subjectIDs)
        rows = T.subjectID == subjectIDs(i);
        thisData = sortrows(T(rows, :), 'subepochIdx');
        if strcmpi(thisData.groupLabel(1), "CP")
            thisColor = cpColor;
        else
            thisColor = controlColor;
        end

        jitter = (thisData.subepochIdx - mean(thisData.subepochIdx)) * 0.045;
        plot(ax, thisData.ageYears + jitter, thisData.alpha_dB, 'o-', ...
            'Color', thisColor, 'LineWidth', 0.9, 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
    end

    hControl = plot(ax, nan, nan, 'o-', 'Color', controlColor, ...
        'LineWidth', 1.2, 'MarkerSize', 5);
    hCP = plot(ax, nan, nan, 'o-', 'Color', cpColor, ...
        'LineWidth', 1.2, 'MarkerSize', 5);

    xlabel(ax, 'Age (years), jittered by sub-epoch');
    ylabel(ax, 'Alpha power (dB)');
    title(ax, 'Sub-epoch alpha estimates per subject');
    legend(ax, [hControl hCP], {'Control', 'CP'}, 'Location', 'best');
    grid(ax, 'on');
    box(ax, 'off');
    set(ax, 'FontName', 'Helvetica', 'FontSize', 11);

    saveas(fig, plotFile);
end
