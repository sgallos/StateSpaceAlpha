function results = SS_age_diff_em(subjectTable, varargin)
%SS_AGE_DIFF_EM Estimate IWP smoothness for the age-difference SSM by EM.
%
% This is the Step 4 version of the StateSpaceAlpha model. It keeps the
% Step 3 observation model:
%
%   x(a) = [f0(a); f0_slope(a); delta(a); delta_slope(a)]
%   y_k  = [1 0 g_k 0] * x(a_k) + noise_k
%
% where g_k = 0 for controls and g_k = 1 for CP, so delta(a) is CP minus
% Control alpha power. The EM loop estimates the two IWP process-noise
% hyperparameters q_f0, q_delta, and sigmaBio. sigmaBio is an additive
% between-subject observation variance:
%
%   V_k = mfdb_var(k) + sigmaBio

    parser = inputParser();
    parser.FunctionName = 'SS_age_diff_em';
    addParameter(parser, 'maxIter', 100, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(parser, 'tolerance', 1e-5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'qInitMultipliers', [0.01 0.1 1 10 100], ...
        @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x > 0));
    addParameter(parser, 'initialProcessNoiseQF0', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'initialProcessNoiseQDelta', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'qFloor', 1e-8, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'qCeiling', 1e6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'initialSigmaBio', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0));
    addParameter(parser, 'sigmaBioFloor', 1e-6, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(parser, 'sigmaBioCeiling', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'initialCovarianceScale', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'minObservationVariance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'useApproximateMstep', false, @(x) islogical(x) || isnumeric(x));
    addParameter(parser, 'verbose', true, @(x) islogical(x) || isnumeric(x));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.useApproximateMstep = logical(opts.useApproximateMstep);
    opts.verbose = logical(opts.verbose);

    data = prepare_model_data(subjectTable, opts);
    qInit = compute_q_init(data.ageYears, data.alpha_dB);

    if isempty(opts.initialProcessNoiseQF0)
        baseInitialQF0 = qInit;
    else
        baseInitialQF0 = opts.initialProcessNoiseQF0;
    end

    if isempty(opts.initialSigmaBio)
        baseInitialSigmaBio = estimate_initial_sigma_bio(data);
    else
        baseInitialSigmaBio = opts.initialSigmaBio;
    end
    baseInitialSigmaBio = clamp_sigma_bio(baseInitialSigmaBio, opts);

    if isempty(opts.initialProcessNoiseQDelta)
        baseInitialQDelta = qInit;
    else
        baseInitialQDelta = opts.initialProcessNoiseQDelta;
    end

    qInitMultipliers = opts.qInitMultipliers(:);
    startResults = cell(numel(qInitMultipliers), 1);
    finalLogLikelihood = nan(numel(qInitMultipliers), 1);
    converged = false(numel(qInitMultipliers), 1);

    if opts.verbose
        fprintf('SS_age_diff_em Step 4 EM setup\n');
        fprintf('Subjects: %d | Control: %d | CP: %d\n', ...
            height(data), sum(data.groupIndicator == 0), sum(data.groupIndicator == 1));
        fprintf('q_init heuristic = %.6g\n', qInit);
        fprintf('initial sigmaBio = %.6g dB^2\n', baseInitialSigmaBio);
        if opts.useApproximateMstep
            fprintf('M-step mode: approximate, without lag-one covariance.\n');
        else
            fprintf('M-step mode: exact, lag-one covariance from joint posterior precision.\n');
        end
    end

    for startIdx = 1:numel(qInitMultipliers)
        multiplier = qInitMultipliers(startIdx);
        startQF0 = clamp_q(baseInitialQF0 * multiplier, opts);
        startQDelta = clamp_q(baseInitialQDelta * multiplier, opts);
        startSigmaBio = baseInitialSigmaBio;

        if opts.verbose
            fprintf('\nEM start %d/%d: q_f0 = %.6g, q_delta = %.6g, sigmaBio = %.6g\n', ...
                startIdx, numel(qInitMultipliers), startQF0, startQDelta, startSigmaBio);
        end

        startResults{startIdx} = run_single_em_start(data, startQF0, startQDelta, startSigmaBio, opts);
        finalLogLikelihood(startIdx) = startResults{startIdx}.finalLogLikelihood;
        converged(startIdx) = startResults{startIdx}.converged;
    end

    [~, bestStartIndex] = max(finalLogLikelihood);
    bestStart = startResults{bestStartIndex};

    results = struct();
    results.modelName = "baseline_delta_4d_iwp_em_step4";
    results.subjectTable = data;
    results.qInit = qInit;
    results.qInitMultipliers = qInitMultipliers;
    results.startResults = startResults;
    results.bestStartIndex = bestStartIndex;
    results.bestStart = bestStart;
    results.trajectory = bestStart.finalFit.trajectory;
    results.q_f0_em = bestStart.q_f0_em;
    results.q_delta_em = bestStart.q_delta_em;
    results.sigmaBio_em = bestStart.sigmaBio_em;
    results.hitQFloor = bestStart.hitQFloor;
    results.hitQCeiling = bestStart.hitQCeiling;
    results.hitSigmaBioFloor = bestStart.hitSigmaBioFloor;
    results.hitSigmaBioCeiling = bestStart.hitSigmaBioCeiling;
    results.finalLogLikelihood = bestStart.finalLogLikelihood;
    results.converged = bestStart.converged;
    results.anyConverged = any(converged);
    results.options = opts;

    if opts.verbose
        fprintf('\nBest EM start: %d\n', bestStartIndex);
        fprintf('q_f0_em = %.6g | q_delta_em = %.6g | sigmaBio_em = %.6g dB^2\n', ...
            results.q_f0_em, results.q_delta_em, results.sigmaBio_em);
        fprintf('Final log-likelihood = %.6f\n', results.finalLogLikelihood);
        fprintf('Converged: %d | iterations: %d\n', ...
            results.converged, results.bestStart.emIterations);
        if results.hitQFloor || results.hitQCeiling || results.hitSigmaBioFloor || results.hitSigmaBioCeiling
            warning('Best EM result hit a hyperparameter boundary. Treat the selected values cautiously.');
        end
    end
end

function startResult = run_single_em_start(data, qF0Start, qDeltaStart, sigmaBioStart, opts)
    maxStoredIter = opts.maxIter + 1;
    qHistory = nan(maxStoredIter, 2);
    sigmaBioHistory = nan(maxStoredIter, 1);
    qNewHistory = nan(opts.maxIter, 2);
    sigmaBioNewHistory = nan(opts.maxIter, 1);
    logLikHistory = nan(maxStoredIter, 1);
    relativeChangeHistory = nan(opts.maxIter, 1);
    mstepDiagnostics = cell(opts.maxIter, 1);

    qF0 = qF0Start;
    qDelta = qDeltaStart;
    sigmaBio = sigmaBioStart;
    converged = false;

    for iter = 1:opts.maxIter
        fit = run_filter_smoother_4d(data, qF0, qDelta, sigmaBio, opts);
        qHistory(iter, :) = [qF0 qDelta];
        sigmaBioHistory(iter) = sigmaBio;
        logLikHistory(iter) = fit.logLikelihood;

        [qF0New, qDeltaNew, sigmaBioNew, thisMstepDiagnostics] = ...
            update_hyperparameters_from_smoother(fit, data, opts);

        qF0New = clamp_q(qF0New, opts);
        qDeltaNew = clamp_q(qDeltaNew, opts);
        sigmaBioNew = clamp_sigma_bio(sigmaBioNew, opts);
        qNewHistory(iter, :) = [qF0New qDeltaNew];
        sigmaBioNewHistory(iter) = sigmaBioNew;
        mstepDiagnostics{iter} = thisMstepDiagnostics;

        denominator = max([qF0, qDelta, sigmaBio], ...
            [opts.qFloor, opts.qFloor, opts.sigmaBioFloor]);
        relativeChange = max(abs([qF0New - qF0, qDeltaNew - qDelta, sigmaBioNew - sigmaBio]) ./ ...
            denominator);
        relativeChangeHistory(iter) = relativeChange;

        if opts.verbose
            fprintf(['  iter %02d: logLik %.6f | q_f0 %.6g -> %.6g | ' ...
                'q_delta %.6g -> %.6g | sigmaBio %.6g -> %.6g | rel %.3g\n'], ...
                iter, fit.logLikelihood, qF0, qF0New, qDelta, qDeltaNew, ...
                sigmaBio, sigmaBioNew, relativeChange);
        end

        qF0 = qF0New;
        qDelta = qDeltaNew;
        sigmaBio = sigmaBioNew;

        if relativeChange < opts.tolerance
            converged = true;
            break;
        end
    end

    finalFit = run_filter_smoother_4d(data, qF0, qDelta, sigmaBio, opts);
    finalIndex = min(iter + 1, maxStoredIter);
    qHistory(finalIndex, :) = [qF0 qDelta];
    sigmaBioHistory(finalIndex) = sigmaBio;
    logLikHistory(finalIndex) = finalFit.logLikelihood;

    validQRows = all(isfinite(qHistory), 2);
    validSigmaRows = isfinite(sigmaBioHistory);
    validLogRows = isfinite(logLikHistory);
    validRelativeRows = isfinite(relativeChangeHistory);

    startResult = struct();
    startResult.initialQF0 = qF0Start;
    startResult.initialQDelta = qDeltaStart;
    startResult.initialSigmaBio = sigmaBioStart;
    startResult.q_f0_em = qF0;
    startResult.q_delta_em = qDelta;
    startResult.sigmaBio_em = sigmaBio;
    startResult.hitQFloor = any(abs([qF0, qDelta] - opts.qFloor) <= eps(opts.qFloor));
    startResult.hitQCeiling = any(abs([qF0, qDelta] - opts.qCeiling) <= eps(opts.qCeiling));
    startResult.hitSigmaBioFloor = abs(sigmaBio - opts.sigmaBioFloor) <= eps(max(opts.sigmaBioFloor, 1));
    startResult.hitSigmaBioCeiling = abs(sigmaBio - opts.sigmaBioCeiling) <= eps(opts.sigmaBioCeiling);
    startResult.finalLogLikelihood = finalFit.logLikelihood;
    startResult.finalFit = finalFit;
    startResult.qHistory = qHistory(validQRows, :);
    startResult.sigmaBioHistory = sigmaBioHistory(validSigmaRows);
    startResult.qNewHistory = qNewHistory(any(isfinite(qNewHistory), 2), :);
    startResult.sigmaBioNewHistory = sigmaBioNewHistory(isfinite(sigmaBioNewHistory));
    startResult.logLikHistory = logLikHistory(validLogRows);
    startResult.relativeChangeHistory = relativeChangeHistory(validRelativeRows);
    startResult.mstepDiagnostics = mstepDiagnostics(~cellfun(@isempty, mstepDiagnostics));
    startResult.converged = converged;
    startResult.emIterations = iter;
    startResult.logLikelihoodIsMonotone = all(diff(startResult.logLikHistory) >= -1e-7);

    if opts.verbose && ~startResult.logLikelihoodIsMonotone
        warning('EM log-likelihood decreased for this start. Check M-step diagnostics.');
    end
end

function fit = run_filter_smoother_4d(data, processNoiseQF0, processNoiseQDelta, sigmaBio, opts)
    ageYears = data.ageYears(:);
    alphaDB = data.alpha_dB(:);
    observationVariance = max(data.mfdb_var(:) + sigmaBio, opts.minObservationVariance);
    groupIndicator = data.groupIndicator(:);

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

    initialState = make_initial_state(alphaDB);
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

    if opts.useApproximateMstep
        lagOneCovariance = nan(stateDimension, stateDimension, nSubjects - 1);
        jointMeanMaxAbsDiff = NaN;
        jointCovMaxAbsDiff = NaN;
    else
        [jointMean, jointCovariance, lagOneCovariance] = ...
            compute_joint_posterior_blocks(data, processNoiseQF0, processNoiseQDelta, ...
            sigmaBio, initialState, initialCovariance, opts);
        jointMarginalCovariance = extract_marginal_covariances(jointCovariance, stateDimension, nSubjects);
        jointMeanMaxAbsDiff = max(abs(jointMean(:) - smoothedState(:)));
        jointCovMaxAbsDiff = max(abs(jointMarginalCovariance(:) - smoothedCovariance(:)));
    end

    trajectory = build_trajectory_table(data, observationVariance, smoothedState, smoothedCovariance);

    fit = struct();
    fit.processNoiseQF0 = processNoiseQF0;
    fit.processNoiseQDelta = processNoiseQDelta;
    fit.sigmaBio = sigmaBio;
    fit.observationVariance = observationVariance;
    fit.trajectory = trajectory;
    fit.predictedState = predictedState;
    fit.filteredState = filteredState;
    fit.smoothedState = smoothedState;
    fit.predictedCovariance = predictedCovariance;
    fit.filteredCovariance = filteredCovariance;
    fit.smoothedCovariance = smoothedCovariance;
    fit.transitionMatrix = transitionMatrix;
    fit.processCovariance = processCovariance;
    fit.observationMatrix = observationMatrix;
    fit.kalmanGain = kalmanGain;
    fit.innovation = innovation;
    fit.innovationVariance = innovationVariance;
    fit.logLikContribution = logLikContribution;
    fit.logLikelihood = sum(logLikContribution);
    fit.lagOneCovariance = lagOneCovariance;
    fit.jointMeanMaxAbsDiff = jointMeanMaxAbsDiff;
    fit.jointCovMaxAbsDiff = jointCovMaxAbsDiff;
end

function [qF0New, qDeltaNew, sigmaBioNew, diagnostics] = update_hyperparameters_from_smoother(fit, data, opts)
    nSubjects = height(data);
    nTransitions = nSubjects - 1;
    contributionF0 = nan(nTransitions, 1);
    contributionDelta = nan(nTransitions, 1);
    rawTraceF0 = nan(nTransitions, 1);
    rawTraceDelta = nan(nTransitions, 1);

    for k = 1:nTransitions
        h = data.ageYears(k + 1) - data.ageYears(k);
        Ak = fit.transitionMatrix(:, :, k);
        baseQ = make_iwp_base_covariance(h);

        Pk = fit.smoothedCovariance(:, :, k);
        Pnext = fit.smoothedCovariance(:, :, k + 1);
        meanInnovation = fit.smoothedState(:, k + 1) - Ak * fit.smoothedState(:, k);

        if opts.useApproximateMstep
            expectedInnovationOuter = Pnext + Ak * Pk * Ak' + ...
                meanInnovation * meanInnovation';
        else
            PkNext = fit.lagOneCovariance(:, :, k);
            expectedInnovationOuter = Pnext + Ak * Pk * Ak' - ...
                PkNext' * Ak' - Ak * PkNext + meanInnovation * meanInnovation';
        end

        expectedInnovationOuter = symmetrize_covariance(expectedInnovationOuter);
        S_f0 = symmetrize_covariance(expectedInnovationOuter(1:2, 1:2));
        S_delta = symmetrize_covariance(expectedInnovationOuter(3:4, 3:4));

        rawTraceF0(k) = trace(baseQ \ S_f0);
        rawTraceDelta(k) = trace(baseQ \ S_delta);
        contributionF0(k) = 0.5 * rawTraceF0(k);
        contributionDelta(k) = 0.5 * rawTraceDelta(k);
    end

    qF0New = sum(max(contributionF0, 0)) / nTransitions;
    qDeltaNew = sum(max(contributionDelta, 0)) / nTransitions;

    residualVarianceTerms = nan(nSubjects, 1);
    for k = 1:nSubjects
        Ck = fit.observationMatrix(k, :);
        smoothedMean = fit.smoothedState(:, k);
        smoothedCovariance = fit.smoothedCovariance(:, :, k);
        innovationSquared = (data.alpha_dB(k) - Ck * smoothedMean)^2;
        smoothedObsVariance = Ck * smoothedCovariance * Ck';
        residualVarianceTerms(k) = innovationSquared + smoothedObsVariance - data.mfdb_var(k);
    end
    sigmaBioRaw = mean(residualVarianceTerms);
    sigmaBioNew = max(sigmaBioRaw, opts.sigmaBioFloor);

    diagnostics = struct();
    diagnostics.contributionF0 = contributionF0;
    diagnostics.contributionDelta = contributionDelta;
    diagnostics.rawTraceF0 = rawTraceF0;
    diagnostics.rawTraceDelta = rawTraceDelta;
    diagnostics.residualVarianceTerms = residualVarianceTerms;
    diagnostics.sigmaBioRaw = sigmaBioRaw;
    diagnostics.anyNegativeContribution = any(contributionF0 < -1e-8) || any(contributionDelta < -1e-8);
end

function [jointMean, jointCovariance, lagOneCovariance] = compute_joint_posterior_blocks( ...
        data, processNoiseQF0, processNoiseQDelta, sigmaBio, initialState, initialCovariance, opts)
    nSubjects = height(data);
    stateDimension = 4;
    totalDimension = stateDimension * nSubjects;
    precisionMatrix = zeros(totalDimension, totalDimension);
    informationVector = zeros(totalDimension, 1);

    idx1 = state_index(1, stateDimension);
    initialPrecision = initialCovariance \ eye(stateDimension);
    precisionMatrix(idx1, idx1) = precisionMatrix(idx1, idx1) + initialPrecision;
    informationVector(idx1) = informationVector(idx1) + initialPrecision * initialState;

    for k = 1:(nSubjects - 1)
        h = data.ageYears(k + 1) - data.ageYears(k);
        Ak = make_iwp_transition_4d(h);
        Qk = make_iwp_process_covariance_4d(h, processNoiseQF0, processNoiseQDelta);
        Qinv = Qk \ eye(stateDimension);

        idxK = state_index(k, stateDimension);
        idxNext = state_index(k + 1, stateDimension);

        precisionMatrix(idxK, idxK) = precisionMatrix(idxK, idxK) + Ak' * Qinv * Ak;
        precisionMatrix(idxK, idxNext) = precisionMatrix(idxK, idxNext) - Ak' * Qinv;
        precisionMatrix(idxNext, idxK) = precisionMatrix(idxNext, idxK) - Qinv * Ak;
        precisionMatrix(idxNext, idxNext) = precisionMatrix(idxNext, idxNext) + Qinv;
    end

    for k = 1:nSubjects
        Ck = [1 0 data.groupIndicator(k) 0];
        observationVariance = max(data.mfdb_var(k) + sigmaBio, opts.minObservationVariance);
        Rinv = 1 / observationVariance;
        idxK = state_index(k, stateDimension);

        precisionMatrix(idxK, idxK) = precisionMatrix(idxK, idxK) + Ck' * Rinv * Ck;
        informationVector(idxK) = informationVector(idxK) + Ck' * Rinv * data.alpha_dB(k);
    end

    precisionMatrix = symmetrize_covariance(precisionMatrix);
    jointMeanVector = precisionMatrix \ informationVector;
    jointCovariance = precisionMatrix \ eye(totalDimension);
    jointCovariance = symmetrize_covariance(jointCovariance);
    jointMean = reshape(jointMeanVector, stateDimension, nSubjects);

    lagOneCovariance = nan(stateDimension, stateDimension, nSubjects - 1);
    for k = 1:(nSubjects - 1)
        idxK = state_index(k, stateDimension);
        idxNext = state_index(k + 1, stateDimension);
        lagOneCovariance(:, :, k) = jointCovariance(idxK, idxNext);
    end
end

function marginalCovariances = extract_marginal_covariances(jointCovariance, stateDimension, nSubjects)
    marginalCovariances = nan(stateDimension, stateDimension, nSubjects);
    for k = 1:nSubjects
        idxK = state_index(k, stateDimension);
        marginalCovariances(:, :, k) = jointCovariance(idxK, idxK);
    end
end

function idx = state_index(k, stateDimension)
    idx = ((k - 1) * stateDimension + 1):(k * stateDimension);
end

function trajectory = build_trajectory_table(data, observationVariance, smoothedState, smoothedCovariance)
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
        data.groupIndicator, ...
        data.ageYears, ...
        data.alpha_dB, ...
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
end

function data = prepare_model_data(subjectTable, opts)
    data = standardize_subject_table(subjectTable);
    validMask = isfinite(data.ageYears) & isfinite(data.alpha_dB) & isfinite(data.mfdb_var);
    data = data(validMask, :);
    data = sortrows(data, {'ageYears', 'subjectID'});

    if height(data) < 4
        error('SS_age_diff_em needs at least 4 valid subjects. Found %d.', height(data));
    end

    data.groupIndicator = make_group_indicator(data.groupLabel);

    ageGaps = diff(data.ageYears);
    if any(ageGaps <= 0)
        error('Subject ages must be strictly increasing after sorting. Duplicate or reversed ages found.');
    end

    nCP = sum(data.groupIndicator == 1);
    nControl = sum(data.groupIndicator == 0);
    if nCP == 0 || nControl == 0
        error('Both CP and Control subjects are required. Found CP=%d, Control=%d.', nCP, nControl);
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

function qInit = compute_q_init(ageYears, alphaDB)
    ageSpan = max(ageYears) - min(ageYears);
    alphaRange = max(alphaDB) - min(alphaDB);
    qInit = (alphaRange^2) / max(ageSpan^3, eps);
end

function initialState = make_initial_state(alphaDB)
    initialState = [median(alphaDB, 'omitnan'); 0; 0; 0];
end

function qValue = clamp_q(qValue, opts)
    qValue = min(max(qValue, opts.qFloor), opts.qCeiling);
end

function sigmaBio = clamp_sigma_bio(sigmaBio, opts)
    sigmaBio = min(max(sigmaBio, opts.sigmaBioFloor), opts.sigmaBioCeiling);
end

function sigmaBio = estimate_initial_sigma_bio(data)
    groupIndicator = make_group_indicator(data.groupLabel);
    ageCentered = data.ageYears - mean(data.ageYears, 'omitnan');
    designMatrix = [ones(height(data), 1), ageCentered, groupIndicator, ageCentered .* groupIndicator];
    beta = designMatrix \ data.alpha_dB;
    residuals = data.alpha_dB - designMatrix * beta;
    sigmaBio = max(var(residuals, 1) - mean(data.mfdb_var, 'omitnan'), 1e-6);
end

function A = make_iwp_transition_4d(ageGapYears)
    block = [1 ageGapYears; 0 1];
    A = blkdiag(block, block);
end

function baseQ = make_iwp_base_covariance(ageGapYears)
    h = ageGapYears;
    baseQ = [h^3 / 3, h^2 / 2; h^2 / 2, h];
end

function Q = make_iwp_process_covariance_4d(ageGapYears, processNoiseQF0, processNoiseQDelta)
    baseQ = make_iwp_base_covariance(ageGapYears);
    Q = blkdiag(processNoiseQF0 * baseQ, processNoiseQDelta * baseQ);
end

function P = symmetrize_covariance(P)
    P = (P + P') / 2;
end
