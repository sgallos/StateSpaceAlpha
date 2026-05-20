function results = SS_age_diff(subjectTable, varargin)
%SS_AGE_DIFF Step 2 pooled age-axis state-space smoother for alpha power.
%
% This Step 2 version deliberately ignores CP/control group labels and fits
% one shared smooth trajectory through all subjects. The scientific
% group-difference model comes later. This file tests the Kalman filter, RTS
% smoother, irregular age gaps, and MFDB observation variances.
%
% Required subjectTable variables, case-insensitive:
%   subjectID   subject identifier
%   groupLabel  CP or Control
%   ageYears    age in years
%   alpha_dB    extracted alpha-band power in dB
%   mfdb_var    per-subject MFDB variance of alpha_dB, in dB^2

    parser = inputParser();
    parser.FunctionName = 'SS_age_diff';
    addParameter(parser, 'qScale', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'processNoiseQ', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'initialCovarianceScale', 100, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'minObservationVariance', 1e-6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'verbose', true, @(x) islogical(x) || isnumeric(x));
    parse(parser, varargin{:});
    opts = parser.Results;

    data = standardize_subject_table(subjectTable);
    data = data(isfinite(data.ageYears) & isfinite(data.alpha_dB) & isfinite(data.mfdb_var), :);
    data = sortrows(data, {'ageYears', 'subjectID'});

    if height(data) < 3
        error('SS_age_diff needs at least 3 valid subjects. Found %d.', height(data));
    end

    ageYears = data.ageYears(:);
    alphaDB = data.alpha_dB(:);
    observationVariance = max(data.mfdb_var(:), opts.minObservationVariance);
    ageGaps = diff(ageYears);

    if any(ageGaps <= 0)
        error('Subject ages must be strictly increasing after sorting. Duplicate or reversed ages found.');
    end

    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange^2) / max(ageSpan^3, eps);

    if isempty(opts.processNoiseQ)
        processNoiseQ = opts.qScale * qInit;
    else
        processNoiseQ = opts.processNoiseQ;
    end

    nSubjects = height(data);
    observationMatrix = [1 0];
    stateDimension = 2;

    predictedState = nan(stateDimension, nSubjects);
    filteredState = nan(stateDimension, nSubjects);
    smoothedState = nan(stateDimension, nSubjects);

    predictedCovariance = nan(stateDimension, stateDimension, nSubjects);
    filteredCovariance = nan(stateDimension, stateDimension, nSubjects);
    smoothedCovariance = nan(stateDimension, stateDimension, nSubjects);
    transitionMatrix = nan(stateDimension, stateDimension, max(nSubjects - 1, 1));
    processCovariance = nan(stateDimension, stateDimension, max(nSubjects - 1, 1));
    kalmanGain = nan(stateDimension, nSubjects);
    innovation = nan(nSubjects, 1);
    innovationVariance = nan(nSubjects, 1);

    initialState = [median(alphaDB, 'omitnan'); 0];
    initialCovariance = opts.initialCovarianceScale * eye(stateDimension);

    for k = 1:nSubjects
        if k == 1
            xPrior = initialState;
            pPrior = initialCovariance;
        else
            h = ageYears(k) - ageYears(k - 1);
            Ak = make_iwp_transition(h);
            Qk = make_iwp_process_covariance(h, processNoiseQ);
            transitionMatrix(:, :, k - 1) = Ak;
            processCovariance(:, :, k - 1) = Qk;

            xPrior = Ak * filteredState(:, k - 1);
            pPrior = Ak * filteredCovariance(:, :, k - 1) * Ak' + Qk;
            pPrior = symmetrize_covariance(pPrior);
        end

        predictedState(:, k) = xPrior;
        predictedCovariance(:, :, k) = pPrior;

        innovation(k) = alphaDB(k) - observationMatrix * xPrior;
        innovationVariance(k) = observationMatrix * pPrior * observationMatrix' + observationVariance(k);

        if innovationVariance(k) <= 0 || ~isfinite(innovationVariance(k))
            error('Invalid innovation variance at subject %d.', k);
        end

        thisKalmanGain = (pPrior * observationMatrix') / innovationVariance(k);
        xPost = xPrior + thisKalmanGain * innovation(k);

        updateMatrix = eye(stateDimension) - thisKalmanGain * observationMatrix;
        pPost = updateMatrix * pPrior * updateMatrix' + ...
            thisKalmanGain * observationVariance(k) * thisKalmanGain';

        filteredState(:, k) = xPost;
        filteredCovariance(:, :, k) = symmetrize_covariance(pPost);
        kalmanGain(:, k) = thisKalmanGain;
    end

    smoothedState(:, nSubjects) = filteredState(:, nSubjects);
    smoothedCovariance(:, :, nSubjects) = filteredCovariance(:, :, nSubjects);

    for k = (nSubjects - 1):-1:1
        Ak = transitionMatrix(:, :, k);
        smootherGain = filteredCovariance(:, :, k) * Ak' / predictedCovariance(:, :, k + 1);

        smoothedState(:, k) = filteredState(:, k) + ...
            smootherGain * (smoothedState(:, k + 1) - predictedState(:, k + 1));
        smoothedCovariance(:, :, k) = filteredCovariance(:, :, k) + ...
            smootherGain * (smoothedCovariance(:, :, k + 1) - predictedCovariance(:, :, k + 1)) * smootherGain';
        smoothedCovariance(:, :, k) = symmetrize_covariance(smoothedCovariance(:, :, k));
    end

    filteredLevelVariance = max(squeeze(filteredCovariance(1, 1, :)), 0);
    smoothedLevelVariance = max(squeeze(smoothedCovariance(1, 1, :)), 0);
    filteredLevelSD = sqrt(filteredLevelVariance);
    smoothedLevelSD = sqrt(smoothedLevelVariance);

    trajectory = table( ...
        data.subjectID, ...
        data.groupLabel, ...
        ageYears, ...
        alphaDB, ...
        observationVariance, ...
        filteredState(1, :).', ...
        filteredLevelSD(:), ...
        smoothedState(1, :).', ...
        smoothedLevelSD(:), ...
        smoothedState(1, :).' - 1.96 * smoothedLevelSD(:), ...
        smoothedState(1, :).' + 1.96 * smoothedLevelSD(:), ...
        filteredState(2, :).', ...
        smoothedState(2, :).', ...
        'VariableNames', {'subjectID','groupLabel','ageYears','alpha_dB','observationVariance', ...
        'filterMean_dB','filterSD_dB','smoothMean_dB','smoothSD_dB', ...
        'smoothCILow_dB','smoothCIHigh_dB','filterSlope_dBPerYear','smoothSlope_dBPerYear'});

    diagnostics = struct();
    diagnostics.qInit = qInit;
    diagnostics.qScale = opts.qScale;
    diagnostics.processNoiseQ = processNoiseQ;
    diagnostics.initialState = initialState;
    diagnostics.initialCovariance = initialCovariance;
    diagnostics.ageGaps = ageGaps;
    diagnostics.innovation = innovation;
    diagnostics.innovationVariance = innovationVariance;
    diagnostics.kalmanGain = kalmanGain;

    results = struct();
    results.modelName = "pooled_2d_iwp_step2";
    results.subjectTable = data;
    results.trajectory = trajectory;
    results.predictedState = predictedState;
    results.filteredState = filteredState;
    results.smoothedState = smoothedState;
    results.predictedCovariance = predictedCovariance;
    results.filteredCovariance = filteredCovariance;
    results.smoothedCovariance = smoothedCovariance;
    results.transitionMatrix = transitionMatrix;
    results.processCovariance = processCovariance;
    results.observationMatrix = observationMatrix;
    results.diagnostics = diagnostics;

    if opts.verbose
        fprintf('SS_age_diff Step 2 pooled IWP fit complete.\n');
        fprintf('Subjects: %d\n', nSubjects);
        fprintf('Age span: %.3f to %.3f years\n', min(ageYears), max(ageYears));
        fprintf('q_init = %.6g | q_scale = %.3g | q = %.6g\n', ...
            qInit, opts.qScale, processNoiseQ);
        fprintf('Observation variance range: %.6g to %.6g dB^2\n', ...
            min(observationVariance), max(observationVariance));
    end
end

function data = standardize_subject_table(subjectTable)
    names = subjectTable.Properties.VariableNames;

    subjectID = get_table_column(subjectTable, names, {'subjectID', 'SubjectID'});
    groupLabel = get_table_column(subjectTable, names, {'groupLabel', 'Group', 'GroupLabel'});
    ageYears = get_table_column(subjectTable, names, {'ageYears', 'AgeYears'});
    alphaDB = get_table_column(subjectTable, names, {'alpha_dB', 'AlphaPower_dB', 'AlphaPowerDB'});
    mfdbVar = get_table_column(subjectTable, names, {'mfdb_var', 'MFDBVar', 'MfdbVar'});

    data = table( ...
        string(subjectID(:)), ...
        string(groupLabel(:)), ...
        double(ageYears(:)), ...
        double(alphaDB(:)), ...
        double(mfdbVar(:)), ...
        'VariableNames', {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var'});
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

function A = make_iwp_transition(ageGapYears)
    A = [1 ageGapYears; 0 1];
end

function Q = make_iwp_process_covariance(ageGapYears, processNoiseQ)
    h = ageGapYears;
    Q = processNoiseQ * [h^3 / 3, h^2 / 2; h^2 / 2, h];
end

function P = symmetrize_covariance(P)
    P = (P + P') / 2;
end
