function results = SS_age_diff_bootstrap_simultaneous(alphaTable, varargin)
%SS_AGE_DIFF_BOOTSTRAP_SIMULTANEOUS Build MFDB bootstrap bands for delta(a).
%
% This Step 5 wrapper can resample MFDB rows only, or do a two-stage
% bootstrap that resamples subjects within group and then draws one MFDB
% alpha-power row for each selected subject. It refits the fixed-q Step 3
% SSM and uses the empirical distribution of
% fitted CP-Control difference trajectories to compute pointwise and
% simultaneous confidence bands.
%
% The smoothness q values are fixed throughout the bootstrap. Do not rerun
% EM inside the bootstrap loop.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    defaultRegistryFile = fullfile(repoRoot, 'references', 'core_inputs', 'mfdb_subject_registry.csv');
    defaultCacheFile = fullfile(repoRoot, 'outputs', 'boot_alpha_dB_cache.mat');

    parser = inputParser();
    parser.FunctionName = 'SS_age_diff_bootstrap_simultaneous';
    addParameter(parser, 'B', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(parser, 'bandRange', [8 13], @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
    addParameter(parser, 'registryFile', defaultRegistryFile, @(x) ischar(x) || isstring(x));
    addParameter(parser, 'cacheFile', defaultCacheFile, @(x) ischar(x) || isstring(x));
    addParameter(parser, 'rebuildCache', false, @(x) islogical(x) || isnumeric(x));
    addParameter(parser, 'rngSeed', 1, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(parser, 'qScaleF0', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'qScaleDelta', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'processNoiseQF0', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'processNoiseQDelta', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
    addParameter(parser, 'biologicalVariance', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
    addParameter(parser, 'bootstrapMode', "two-stage", @(x) ischar(x) || isstring(x));
    addParameter(parser, 'duplicateAgeJitter', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'progressEvery', 100, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(parser, 'verbose', true, @(x) islogical(x) || isnumeric(x));
    parse(parser, varargin{:});
    opts = parser.Results;
    opts.registryFile = string(opts.registryFile);
    opts.cacheFile = string(opts.cacheFile);
    opts.rebuildCache = logical(opts.rebuildCache);
    opts.verbose = logical(opts.verbose);
    opts.B = round(opts.B);
    opts.progressEvery = round(opts.progressEvery);
    opts.bandRange = double(opts.bandRange(:).');
    opts.bootstrapMode = lower(string(opts.bootstrapMode));
    if ~ismember(opts.bootstrapMode, ["two-stage", "mfdb-only"])
        error('bootstrapMode must be "two-stage" or "mfdb-only".');
    end

    data = prepare_model_data(alphaTable);
    qInit = compute_q_init(data.ageYears, data.alpha_dB);

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

    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    cache = load_or_build_boot_alpha_cache(data, opts);
    validate_cache_alignment(data, cache, opts);

    nSubjects = height(data);
    nAvailableBoot = size(cache.bootAlphaDB, 2);

    originalFit = SS_age_diff(data, ...
        'processNoiseQF0', processNoiseQF0, ...
        'processNoiseQDelta', processNoiseQDelta, ...
        'biologicalVariance', opts.biologicalVariance, ...
        'verbose', false);
    ageGrid = originalFit.trajectory.ageYears(:).';
    originalDelta = originalFit.trajectory.deltaMean_dB(:).';

    deltaBootstrap = nan(opts.B, numel(ageGrid));
    subjectResampleIndices = nan(opts.B, nSubjects);
    mfdbResampleIndices = nan(opts.B, nSubjects);
    failedBootstrap = false(opts.B, 1);
    failureMessages = strings(opts.B, 1);

    if opts.verbose
        fprintf('Step 5 bootstrap: B = %d, subjects = %d, MFDB rows per subject = %d\n', ...
            opts.B, nSubjects, nAvailableBoot);
        fprintf('Bootstrap mode: %s\n', opts.bootstrapMode);
        fprintf('Using fixed q_f0 = %.6g, q_delta = %.6g, biological variance = %.6g dB^2\n', ...
            processNoiseQF0, processNoiseQDelta, opts.biologicalVariance);
    end

    controlIdx = find(data.groupIndicator == 0);
    cpIdx = find(data.groupIndicator == 1);

    for b = 1:opts.B
        switch opts.bootstrapMode
            case "two-stage"
                bootControlIdx = controlIdx(randi(numel(controlIdx), [numel(controlIdx), 1]));
                bootCpIdx = cpIdx(randi(numel(cpIdx), [numel(cpIdx), 1]));
                bootSubjectIdx = [bootControlIdx; bootCpIdx];

            case "mfdb-only"
                bootSubjectIdx = (1:nSubjects).';
        end

        thisDraw = randi(nAvailableBoot, [numel(bootSubjectIdx), 1]);
        subjectResampleIndices(b, 1:numel(bootSubjectIdx)) = bootSubjectIdx;
        mfdbResampleIndices(b, 1:numel(bootSubjectIdx)) = thisDraw;

        resampledTable = data(bootSubjectIdx, {'subjectID','groupLabel','ageYears','alpha_dB','mfdb_var'});
        linearIndex = sub2ind(size(cache.bootAlphaDB), bootSubjectIdx, thisDraw);
        resampledTable.alpha_dB = cache.bootAlphaDB(linearIndex);
        if opts.bootstrapMode == "two-stage"
            resampledTable = jitter_duplicate_ages(resampledTable, opts.duplicateAgeJitter);
        end

        try
            bootFit = SS_age_diff(resampledTable, ...
                'processNoiseQF0', processNoiseQF0, ...
                'processNoiseQDelta', processNoiseQDelta, ...
                'biologicalVariance', opts.biologicalVariance, ...
                'verbose', false);
            deltaBootstrap(b, :) = interp1(bootFit.trajectory.ageYears, ...
                bootFit.trajectory.deltaMean_dB, ageGrid, 'linear', 'extrap');
        catch ME
            failedBootstrap(b) = true;
            failureMessages(b) = string(ME.message);
        end

        if opts.verbose && (mod(b, opts.progressEvery) == 0 || b == opts.B)
            fprintf('Bootstrap iteration %d / %d\n', b, opts.B);
        end
    end

    validBootstrapMask = ~failedBootstrap & all(isfinite(deltaBootstrap), 2);
    if ~any(validBootstrapMask)
        error('All bootstrap SSM fits failed.');
    end

    deltaBootstrapValid = deltaBootstrap(validBootstrapMask, :);
    deltaMean = mean(deltaBootstrapValid, 1);
    deltaSD = std(deltaBootstrapValid, 0, 1);
    deltaSDSafe = max(deltaSD, eps);

    deltaPointwiseLow = deltaMean - 1.96 * deltaSD;
    deltaPointwiseHigh = deltaMean + 1.96 * deltaSD;

    standardizedDeviation = abs(deltaBootstrapValid - deltaMean) ./ deltaSDSafe;
    tBootstrap = max(standardizedDeviation, [], 2);
    tStar = prctile(tBootstrap, 95);

    deltaSimultaneousLow = deltaMean - tStar * deltaSD;
    deltaSimultaneousHigh = deltaMean + tStar * deltaSD;

    pointwiseSigMask = deltaPointwiseLow > 0 | deltaPointwiseHigh < 0;
    simultaneousSigMask = deltaSimultaneousLow > 0 | deltaSimultaneousHigh < 0;
    bootstrapBias = deltaMean - originalDelta;

    bandSummary = table( ...
        ageGrid(:), ...
        originalDelta(:), ...
        deltaMean(:), ...
        deltaSD(:), ...
        bootstrapBias(:), ...
        deltaPointwiseLow(:), ...
        deltaPointwiseHigh(:), ...
        deltaSimultaneousLow(:), ...
        deltaSimultaneousHigh(:), ...
        pointwiseSigMask(:), ...
        simultaneousSigMask(:), ...
        'VariableNames', {'ageYears','originalDelta_dB','bootstrapMeanDelta_dB', ...
        'bootstrapSD_dB','bootstrapBias_dB','pointwiseLow_dB','pointwiseHigh_dB', ...
        'simultaneousLow_dB','simultaneousHigh_dB','pointwiseExcludesZero', ...
        'simultaneousExcludesZero'});

    diagnostics = struct();
    diagnostics.qInit = qInit;
    diagnostics.processNoiseQF0 = processNoiseQF0;
    diagnostics.processNoiseQDelta = processNoiseQDelta;
    diagnostics.biologicalVariance = opts.biologicalVariance;
    diagnostics.bootstrapMode = opts.bootstrapMode;
    diagnostics.nRequestedBootstraps = opts.B;
    diagnostics.nValidBootstraps = sum(validBootstrapMask);
    diagnostics.nFailedBootstraps = sum(failedBootstrap);
    diagnostics.failureMessages = failureMessages(failedBootstrap);
    diagnostics.tStar = tStar;
    diagnostics.tBootstrap = tBootstrap;
    diagnostics.maxAbsBootstrapBias = max(abs(bootstrapBias));
    diagnostics.pointwiseAnyExcludesZero = any(pointwiseSigMask);
    diagnostics.simultaneousAnyExcludesZero = any(simultaneousSigMask);
    diagnostics.cacheFile = opts.cacheFile;
    diagnostics.registryFile = opts.registryFile;
    diagnostics.bandRange = opts.bandRange;
    diagnostics.rngSeed = opts.rngSeed;

    results = struct();
    results.modelName = "baseline_delta_4d_iwp_mfdb_simultaneous_bootstrap_step5";
    results.subjectTable = data;
    results.originalFit = originalFit;
    results.cache = cache;
    results.ageGrid = ageGrid;
    results.deltaBootstrap = deltaBootstrap;
    results.validBootstrapMask = validBootstrapMask;
    results.subjectResampleIndices = subjectResampleIndices;
    results.mfdbResampleIndices = mfdbResampleIndices;
    results.deltaMean = deltaMean;
    results.deltaSD = deltaSD;
    results.deltaPointwiseLow = deltaPointwiseLow;
    results.deltaPointwiseHigh = deltaPointwiseHigh;
    results.deltaSimultaneousLow = deltaSimultaneousLow;
    results.deltaSimultaneousHigh = deltaSimultaneousHigh;
    results.pointwiseSigMask = pointwiseSigMask;
    results.simultaneousSigMask = simultaneousSigMask;
    results.bandSummary = bandSummary;
    results.diagnostics = diagnostics;

    if opts.verbose
        fprintf('Bootstrap valid fits: %d / %d\n', diagnostics.nValidBootstraps, opts.B);
        fprintf('Simultaneous critical value T_star = %.3f\n', tStar);
        fprintf('Max abs bootstrap bias vs original fit = %.3f dB\n', diagnostics.maxAbsBootstrapBias);
        fprintf('Pointwise excludes zero anywhere: %d\n', diagnostics.pointwiseAnyExcludesZero);
        fprintf('Simultaneous excludes zero anywhere: %d\n', diagnostics.simultaneousAnyExcludesZero);
    end
end

function cache = load_or_build_boot_alpha_cache(data, opts)
    cacheFile = opts.cacheFile;
    if isfile(cacheFile) && ~opts.rebuildCache
        loaded = load(cacheFile);
        cache = loaded.cache;
        return;
    end

    if opts.verbose
        fprintf('Building MFDB alpha bootstrap cache from registry:\n  %s\n', opts.registryFile);
    end

    if ~isfile(opts.registryFile)
        error('Registry file not found: %s', opts.registryFile);
    end

    registry = readtable(opts.registryFile, 'TextType', 'string');
    nSubjects = height(data);
    bootAlphaDB = [];
    filePaths = strings(nSubjects, 1);
    nBootRows = nan(nSubjects, 1);

    for k = 1:nSubjects
        matchIdx = find(strcmp(registry.subjectID, data.subjectID(k)), 1, 'first');
        if isempty(matchIdx)
            error('Subject %s was not found in registry.', data.subjectID(k));
        end

        thisFile = registry.filePath(matchIdx);
        filePaths(k) = thisFile;
        if ~isfile(thisFile)
            error('MFDB file not found for %s:\n%s', data.subjectID(k), thisFile);
        end

        S = load(thisFile, 'freq', 'S_subject_boot_dB');
        freq = double(S.freq(:).');
        subjectBootDB = double(S.S_subject_boot_dB);
        if size(subjectBootDB, 2) ~= numel(freq) && size(subjectBootDB, 1) == numel(freq)
            subjectBootDB = subjectBootDB.';
        end

        bandMask = freq >= opts.bandRange(1) & freq <= opts.bandRange(2);
        if ~any(bandMask)
            error('No frequency bins found in %.2f-%.2f Hz for %s.', ...
                opts.bandRange(1), opts.bandRange(2), data.subjectID(k));
        end

        bootLinear = trapz(freq(bandMask), 10.^(subjectBootDB(:, bandMask) / 10), 2);
        thisBootAlphaDB = 10 * log10(max(bootLinear, eps));

        if isempty(bootAlphaDB)
            bootAlphaDB = nan(nSubjects, numel(thisBootAlphaDB));
        elseif size(bootAlphaDB, 2) ~= numel(thisBootAlphaDB)
            error('Subject %s has %d bootstrap rows; expected %d.', ...
                data.subjectID(k), numel(thisBootAlphaDB), size(bootAlphaDB, 2));
        end

        bootAlphaDB(k, :) = thisBootAlphaDB(:).';
        nBootRows(k) = numel(thisBootAlphaDB);

        if opts.verbose
            fprintf('  cached %s (%d MFDB rows)\n', data.subjectID(k), nBootRows(k));
        end
    end

    cache = struct();
    cache.bootAlphaDB = bootAlphaDB;
    cache.subjectID = data.subjectID;
    cache.groupLabel = data.groupLabel;
    cache.groupIndicator = data.groupIndicator;
    cache.ageYears = data.ageYears;
    cache.filePath = filePaths;
    cache.nBootRows = nBootRows;
    cache.bandRange = opts.bandRange;
    cache.registryFile = opts.registryFile;
    cache.createdOn = string(datetime('now'));

    cacheDir = fileparts(cacheFile);
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end
    save(cacheFile, 'cache');
end

function validate_cache_alignment(data, cache, opts)
    if ~isequal(string(cache.subjectID(:)), string(data.subjectID(:)))
        error('Cache subject order does not match the alpha table. Rebuild the cache.');
    end
    if any(abs(cache.ageYears(:) - data.ageYears(:)) > 1e-8)
        error('Cache ages do not match the alpha table. Rebuild the cache.');
    end
    if ~isequal(string(cache.groupLabel(:)), string(data.groupLabel(:)))
        error('Cache group labels do not match the alpha table. Rebuild the cache.');
    end
    if ~isequal(double(cache.bandRange(:).'), double(opts.bandRange(:).'))
        error('Cache band range does not match requested band range. Rebuild the cache.');
    end
end

function T = jitter_duplicate_ages(T, jitterAmount)
    [~, sortIdx] = sortrows(T, {'ageYears', 'subjectID'});
    T = T(sortIdx, :);
    roundedAges = round(T.ageYears, 10);
    uniqueAges = unique(roundedAges, 'stable');

    for i = 1:numel(uniqueAges)
        duplicateIdx = find(roundedAges == uniqueAges(i));
        if numel(duplicateIdx) > 1
            offsets = linspace(-jitterAmount, jitterAmount, numel(duplicateIdx)).';
            T.ageYears(duplicateIdx) = T.ageYears(duplicateIdx) + offsets;
        end
    end
end

function data = prepare_model_data(subjectTable)
    data = standardize_subject_table(subjectTable);
    validMask = isfinite(data.ageYears) & isfinite(data.alpha_dB) & isfinite(data.mfdb_var);
    data = data(validMask, :);
    data = sortrows(data, {'ageYears', 'subjectID'});
    data.groupIndicator = make_group_indicator(data.groupLabel);

    if any(diff(data.ageYears) <= 0)
        error('Subject ages must be strictly increasing after sorting.');
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
