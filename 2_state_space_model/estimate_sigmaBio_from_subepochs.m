function [sigmaBioFixed, diagnostics] = estimate_sigmaBio_from_subepochs(subepochInput, varargin)
%ESTIMATE_SIGMABIO_FROM_SUBEPOCHS Estimate fixed between-subject variance.
%
% This helper uses the 76-row 4 x 30 s sub-epoch alpha table to estimate a
% fixed subject-level variance component for the SSM. It first collapses each
% subject to a mean alpha value, then estimates the between-subject variance
% under one of several explicit estimator choices.
%
% The returned sigmaBioFixed is intended to be held fixed while EM estimates
% only q_f0 and q_delta.

    parser = inputParser();
    parser.FunctionName = 'estimate_sigmaBio_from_subepochs';
    addParameter(parser, 'Estimator', "trend_residual_sem", @(x) isstring(x) || ischar(x));
    addParameter(parser, 'TrendModel', "linear", @(x) isstring(x) || ischar(x));
    addParameter(parser, 'SigmaBioFloor', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(parser, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.Estimator = lower(string(opts.Estimator));
    opts.TrendModel = lower(string(opts.TrendModel));
    opts.Verbose = logical(opts.Verbose);

    if istable(subepochInput)
        subepochTable = subepochInput;
        sourceFile = "";
    else
        sourceFile = string(subepochInput);
        subepochTable = readtable(sourceFile, 'TextType', 'string');
    end

    subepochTable = standardize_subepoch_table(subepochTable);
    subjectTable = build_subject_summary(subepochTable);
    [varianceSource, varianceDescription, residuals] = ...
        make_sigma_bio_variance_source(subjectTable, opts);

    validVarianceSource = varianceSource(isfinite(varianceSource));
    if numel(validVarianceSource) < 4
        error('Need at least 4 valid subject values to estimate sigmaBio.');
    end

    empiricalVarianceOfSubjectMeans = var(validVarianceSource, 0);
    measurementVarianceToSubtract = choose_measurement_variance(subjectTable, opts.Estimator);
    sigmaBioRaw = empiricalVarianceOfSubjectMeans - measurementVarianceToSubtract;
    sigmaBioFixed = max(sigmaBioRaw, opts.SigmaBioFloor);

    subjectTable.trendResidual_dB = residuals;

    diagnostics = struct();
    diagnostics.sourceFile = sourceFile;
    diagnostics.estimator = opts.Estimator;
    diagnostics.trendModel = opts.TrendModel;
    diagnostics.varianceDescription = varianceDescription;
    diagnostics.subjectTable = subjectTable;
    diagnostics.empiricalVarianceOfSubjectMeans = empiricalVarianceOfSubjectMeans;
    diagnostics.meanMeasurementVarianceOfSubjectMean = measurementVarianceToSubtract;
    diagnostics.sigmaBioRaw = sigmaBioRaw;
    diagnostics.sigmaBioFixed = sigmaBioFixed;
    diagnostics.sigmaBioFloor = opts.SigmaBioFloor;

    if opts.Verbose
        fprintf('Fixed sigmaBio estimate from sub-epochs\n');
        fprintf('  estimator: %s\n', opts.Estimator);
        fprintf('  trend model: %s\n', opts.TrendModel);
        fprintf('  variance source: %s\n', varianceDescription);
        fprintf('  empirical variance of subject source: %.6f dB^2\n', ...
            empiricalVarianceOfSubjectMeans);
        fprintf('  variance subtracted: %.6f dB^2\n', measurementVarianceToSubtract);
        fprintf('  raw sigmaBio: %.6f dB^2\n', sigmaBioRaw);
        fprintf('  fixed sigmaBio: %.6f dB^2\n', sigmaBioFixed);
    end
end

function [varianceSource, varianceDescription, residuals] = ...
        make_sigma_bio_variance_source(subjectTable, opts)
    residuals = remove_group_age_trend(subjectTable, opts.TrendModel);

    switch opts.Estimator
        case {"trend_residual_mfdb", "trend_residual_sem"}
            varianceSource = residuals;
            varianceDescription = "group-age trend residuals";

        case "between_only"
            varianceSource = subjectTable.subjectMeanAlpha_dB;
            varianceDescription = "raw subject means, no group or age trend removed";

        otherwise
            error(['Estimator must be "trend_residual_sem", "trend_residual_mfdb", ' ...
                'or "between_only".']);
    end
end

function measurementVariance = choose_measurement_variance(subjectTable, estimator)
    switch estimator
        case "trend_residual_mfdb"
            measurementVariance = mean(subjectTable.measurementVarianceOfMean, 'omitnan');

        case {"trend_residual_sem", "between_only"}
            measurementVariance = mean(subjectTable.empiricalSemVariance, 'omitnan');

        otherwise
            error('Unknown sigmaBio estimator "%s".', estimator);
    end
end

function T = standardize_subepoch_table(T)
    names = T.Properties.VariableNames;
    subjectID = get_table_column(T, names, {'subjectID', 'SubjectID'});
    subepochIdx = get_table_column(T, names, {'subepochIdx', 'SubEpochIdx', 'subEpochIdx'});
    groupLabel = get_table_column(T, names, {'groupLabel', 'Group', 'GroupLabel'});
    ageYears = get_table_column(T, names, {'ageYears', 'AgeYears'});
    alphaDB = get_table_column(T, names, {'alpha_dB', 'AlphaPower_dB', 'AlphaPowerDB'});
    mfdbVar = get_table_column(T, names, {'mfdb_var', 'MFDBVar', 'MfdbVar'});

    T = table(string(subjectID(:)), double(subepochIdx(:)), string(groupLabel(:)), ...
        double(ageYears(:)), double(alphaDB(:)), double(mfdbVar(:)), ...
        'VariableNames', {'subjectID','subepochIdx','groupLabel','ageYears','alpha_dB','mfdb_var'});
    T = sortrows(T, {'ageYears','subjectID','subepochIdx'});
end

function subjectTable = build_subject_summary(subepochTable)
    [G, subjectIDs] = findgroups(subepochTable.subjectID);
    nSubjects = numel(subjectIDs);

    subjectTable = table('Size', [nSubjects, 10], ...
        'VariableTypes', {'string','string','double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'subjectID','groupLabel','ageYears','nSubEpochs','subjectMeanAlpha_dB', ...
        'subjectAlphaSD_dB','subjectAlphaRange_dB','meanSubepochMfdbVariance', ...
        'measurementVarianceOfMean','empiricalSemVariance'});

    for i = 1:nSubjects
        rows = G == i;
        thisAlpha = subepochTable.alpha_dB(rows);
        thisVar = subepochTable.mfdb_var(rows);
        nRows = sum(rows);

        subjectTable.subjectID(i) = subjectIDs(i);
        subjectTable.groupLabel(i) = subepochTable.groupLabel(find(rows, 1));
        subjectTable.ageYears(i) = subepochTable.ageYears(find(rows, 1));
        subjectTable.nSubEpochs(i) = nRows;
        validAlpha = thisAlpha(isfinite(thisAlpha));
        subjectTable.subjectMeanAlpha_dB(i) = mean(validAlpha);
        subjectTable.subjectAlphaSD_dB(i) = std(validAlpha, 0);
        subjectTable.subjectAlphaRange_dB(i) = max(thisAlpha) - min(thisAlpha);
        subjectTable.meanSubepochMfdbVariance(i) = mean(thisVar, 'omitnan');
        subjectTable.measurementVarianceOfMean(i) = sum(thisVar, 'omitnan') / (nRows ^ 2);
        if numel(validAlpha) > 1
            subjectTable.empiricalSemVariance(i) = var(validAlpha, 0) / numel(validAlpha);
        else
            subjectTable.empiricalSemVariance(i) = 0;
        end
    end
end

function residuals = remove_group_age_trend(subjectTable, trendModel)
    residuals = nan(height(subjectTable), 1);
    groupLabels = unique(subjectTable.groupLabel, 'stable');
    degree = trend_model_degree(trendModel);

    for i = 1:numel(groupLabels)
        rows = subjectTable.groupLabel == groupLabels(i);
        age = subjectTable.ageYears(rows);
        y = subjectTable.subjectMeanAlpha_dB(rows);
        thisDegree = min(degree, numel(age) - 1);
        if thisDegree < 0
            continue;
        end

        ageCenter = mean(age, 'omitnan');
        ageCentered = age - ageCenter;
        coeff = polyfit(ageCentered, y, thisDegree);
        fitted = polyval(coeff, ageCentered);
        residuals(rows) = y - fitted;
    end
end

function degree = trend_model_degree(trendModel)
    switch trendModel
        case "constant"
            degree = 0;
        case "linear"
            degree = 1;
        case "quadratic"
            degree = 2;
        otherwise
            error('TrendModel must be "constant", "linear", or "quadratic".');
    end
end

function values = get_table_column(T, allNames, candidates)
    values = [];
    for i = 1:numel(candidates)
        idx = find(strcmp(allNames, candidates{i}), 1, 'first');
        if ~isempty(idx)
            values = T.(allNames{idx});
            return;
        end
    end
    error('Missing required table variable. Tried: %s', strjoin(candidates, ', '));
end
