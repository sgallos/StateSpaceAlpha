function results = SS_age_diff(subjectTable, varargin)
%SS_AGE_DIFF Step 3 4-D age-varying CP-control difference model.
%
% This model estimates a baseline alpha-power age trajectory and a
% CP-control difference trajectory:
%
%   x(a) = [f0(a); f0_slope(a); delta(a); delta_slope(a)]
%
% where delta(a) is CP minus Control. Group labels enter through the
% per-subject observation matrix C_k = [1 0 g_k 0], with g_k = 0 for
% Control and g_k = 1 for CP.
%
% Required subjectTable variables, case-insensitive:
%   subjectID   subject identifier
%   groupLabel  CP or Control
%   ageYears    age in years
%   alpha_dB    extracted alpha-band power in dB
%   mfdb_var    per-subject MFDB variance of alpha_dB, in dB^2

    parser = inputParser();
    parser.FunctionName = 'SS_age_diff';
    addParameter(parser, 'qScaleF0', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'qScaleDelta', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'processNoiseQF0', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'processNoiseQDelta', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'initialCovarianceScale', 100, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'minObservationVariance', 1e-6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'biologicalVariance', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(parser, 'verbose', true, @(x) islogical(x) || isnumeric(x));
    parse(parser, varargin{:});
    opts = parser.Results;

    data = standardize_subject_table(subjectTable);
    data = data(isfinite(data.ageYears) & isfinite(data.alpha_dB) & isfinite(data.mfdb_var), :);
    data = sortrows(data, {'ageYears', 'subjectID'});

    if height(data) < 4
        error('SS_age_diff needs at least 4 valid subjects. Found %d.', height(data));
    end

    ageYears = data.ageYears(:);
    alphaDB = data.alpha_dB(:);
    mfdbObservationVariance = data.mfdb_var(:);
    biologicalVariance = opts.biologicalVariance;
    observationVariance = max(mfdbObservationVariance + biologicalVariance, opts.minObservationVariance);
    groupIndicator = make_group_indicator(data.groupLabel);
    ageGaps = diff(ageYears);

    if any(ageGaps <= 0)
        error('Subject ages must be strictly increasing after sorting. Duplicate or reversed ages found.');
    end

    nCP = sum(groupIndicator == 1);
    nControl = sum(groupIndicator == 0);
    if nCP == 0 || nControl == 0
        error('Both CP and Control subjects are required. Found CP=%d, Control=%d.', nCP, nControl);
    end

    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange^2) / max(ageSpan^3, eps);

    if isempty(opts.processNoiseQF0)
        processNoiseQF0 = opts.qScaleF0 * qInit;
    else
        processNoiseQF0 = opts.processNoiseQF0;
    end

    if isempty(opts.processNoiseQDelta)
        processNoiseQDelta = opts.qScaleDelta * qInit;
    else
        processNoiseQDelta = opts.processNoiseQDelta;
    end

    nSubjects = height(data);
    stateDimension = 4;

    predictedState = nan(stateDimension, nSubjects);
    filteredState = nan(stateDimension, nSubjects);
    smoothedState = nan(stateDimension, nSubjects);

    predictedCovariance = nan(stateDimension, stateDimension, nSubjects);
    filteredCovariance = nan(stateDimension, stateDimension, nSubjects);
    smoothedCovariance = nan(stateDimension, stateDimension, nSubjects);
    transitionMatrix = nan(stateDimension, stateDimension, max(nSubjects - 1, 1));
    processCovariance = nan(stateDimension, stateDimension, max(nSubjects - 1, 1));
    observationMatrix = nan(nSubjects, stateDimension);
    kalmanGain = nan(stateDimension, nSubjects);
    innovation = nan(nSubjects, 1);
    innovationVariance = nan(nSubjects, 1);
    logLikContribution = nan(nSubjects, 1);

    initialState = [median(alphaDB, 'omitnan'); 0; 0; 0];
    initialCovariance = opts.initialCovarianceScale * eye(stateDimension);

    for k = 1:nSubjects
        if k == 1
            xPrior = initialState;
            pPrior = initialCovariance;
        else
            h = ageYears(k) - ageYears(k - 1);
            Ak = make_iwp_transition_4d(h);
            Qk = make_iwp_process_covariance_4d(h, processNoiseQF0, processNoiseQDelta);
            transitionMatrix(:, :, k - 1) = Ak;
            processCovariance(:, :, k - 1) = Qk;

            xPrior = Ak * filteredState(:, k - 1);
            pPrior = Ak * filteredCovariance(:, :, k - 1) * Ak' + Qk;
            pPrior = symmetrize_covariance(pPrior);
        end

        Ck = [1 0 groupIndicator(k) 0];
        observationMatrix(k, :) = Ck;

        predictedState(:, k) = xPrior;
        predictedCovariance(:, :, k) = pPrior;

        innovation(k) = alphaDB(k) - Ck * xPrior;
        innovationVariance(k) = Ck * pPrior * Ck' + observationVariance(k);

        if innovationVariance(k) <= 0 || ~isfinite(innovationVariance(k))
            error('Invalid innovation variance at subject %d.', k);
        end
        logLikContribution(k) = -0.5 * (log(2 * pi * innovationVariance(k)) + ...
            innovation(k)^2 / innovationVariance(k));

        thisKalmanGain = (pPrior * Ck') / innovationVariance(k);
        xPost = xPrior + thisKalmanGain * innovation(k);

        updateMatrix = eye(stateDimension) - thisKalmanGain * Ck;
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

    baselineMean = smoothedState(1, :).';
    deltaMean = smoothedState(3, :).';
    cpMean = baselineMean + deltaMean;

    baselineVariance = max(squeeze(smoothedCovariance(1, 1, :)), 0);
    deltaVariance = max(squeeze(smoothedCovariance(3, 3, :)), 0);
    cpVariance = max(squeeze(smoothedCovariance(1, 1, :) + ...
        2 * smoothedCovariance(1, 3, :) + smoothedCovariance(3, 3, :)), 0);

    baselineSD = sqrt(baselineVariance);
    deltaSD = sqrt(deltaVariance);
    cpSD = sqrt(cpVariance);

    trajectory = table( ...
        data.subjectID, ...
        data.groupLabel, ...
        groupIndicator, ...
        ageYears, ...
        alphaDB, ...
        observationVariance, ...
        baselineMean, ...
        baselineSD(:), ...
        baselineMean - 1.96 * baselineSD(:), ...
        baselineMean + 1.96 * baselineSD(:), ...
        deltaMean, ...
        deltaSD(:), ...
        deltaMean - 1.96 * deltaSD(:), ...
        deltaMean + 1.96 * deltaSD(:), ...
        cpMean, ...
        cpSD(:), ...
        cpMean - 1.96 * cpSD(:), ...
        cpMean + 1.96 * cpSD(:), ...
        smoothedState(2, :).', ...
        smoothedState(4, :).', ...
        'VariableNames', {'subjectID','groupLabel','groupIndicator','ageYears','alpha_dB', ...
        'observationVariance','baselineMean_dB','baselineSD_dB','baselineCILow_dB', ...
        'baselineCIHigh_dB','deltaMean_dB','deltaSD_dB','deltaCILow_dB', ...
        'deltaCIHigh_dB','cpMean_dB','cpSD_dB','cpCILow_dB','cpCIHigh_dB', ...
        'baselineSlope_dBPerYear','deltaSlope_dBPerYear'});

    diagnostics = struct();
    diagnostics.qInit = qInit;
    diagnostics.qScaleF0 = opts.qScaleF0;
    diagnostics.qScaleDelta = opts.qScaleDelta;
    diagnostics.processNoiseQF0 = processNoiseQF0;
    diagnostics.processNoiseQDelta = processNoiseQDelta;
    diagnostics.biologicalVariance = biologicalVariance;
    diagnostics.mfdbObservationVariance = mfdbObservationVariance;
    diagnostics.totalObservationVariance = observationVariance;
    diagnostics.initialState = initialState;
    diagnostics.initialCovariance = initialCovariance;
    diagnostics.ageGaps = ageGaps;
    diagnostics.groupIndicator = groupIndicator;
    diagnostics.nCP = nCP;
    diagnostics.nControl = nControl;
    diagnostics.innovation = innovation;
    diagnostics.innovationVariance = innovationVariance;
    diagnostics.logLikContribution = logLikContribution;
    diagnostics.logLikelihood = sum(logLikContribution);
    diagnostics.kalmanGain = kalmanGain;

    results = struct();
    results.modelName = "baseline_delta_4d_iwp_step3";
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
    results.logLikelihood = diagnostics.logLikelihood;

    if opts.verbose
        fprintf('SS_age_diff Step 3 baseline/delta IWP fit complete.\n');
        fprintf('Subjects: %d | Control: %d | CP: %d\n', nSubjects, nControl, nCP);
        fprintf('Age span: %.3f to %.3f years\n', min(ageYears), max(ageYears));
        fprintf('q_init = %.6g | q_f0 = %.6g | q_delta = %.6g\n', ...
            qInit, processNoiseQF0, processNoiseQDelta);
        fprintf('Biological variance = %.6g dB^2\n', biologicalVariance);
        fprintf('Total observation variance range: %.6g to %.6g dB^2\n', ...
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

function groupIndicator = make_group_indicator(groupLabel)
    groupLabel = string(groupLabel(:));
    isCP = strcmpi(groupLabel, "CP") | strcmpi(groupLabel, "Patient");
    isControl = strcmpi(groupLabel, "Control") | strcmpi(groupLabel, "CN");

    unknownMask = ~(isCP | isControl);
    if any(unknownMask)
        error('Unknown group labels: %s', strjoin(unique(groupLabel(unknownMask)), ', '));
    end

    groupIndicator = double(isCP);
end

function A = make_iwp_transition_4d(ageGapYears)
    block = [1 ageGapYears; 0 1];
    A = blkdiag(block, block);
end

function Q = make_iwp_process_covariance_4d(ageGapYears, processNoiseQF0, processNoiseQDelta)
    h = ageGapYears;
    block = [h^3 / 3, h^2 / 2; h^2 / 2, h];
    Q = blkdiag(processNoiseQF0 * block, processNoiseQDelta * block);
end

function P = symmetrize_covariance(P)
    P = (P + P') / 2;
end
