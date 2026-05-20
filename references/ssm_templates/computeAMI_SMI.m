function [AMI, SMI] = computeAMI_SMI(data, fs, threshold_TIME_L, threshold_TIME_H, varargin)
% Function to compute Alpha Modulation Index (AMI) and Slow Modulation Index (SMI)
% Based on Adam et al. 2023 - Modulatory dynamics mark the transition between 
% anesthetic states of unconsciousness
%
% Inputs:
%   data - raw EEG signal (column vector)
%   fs - sampling frequency (Hz)
%   threshold_TIME_L - start time for calibration window (seconds)
%   threshold_TIME_H - end time for calibration window (seconds)
%
% Optional inputs (name-value pairs):
%   'AlphaThresholdType' - 'percentile' (default) or 'absolute'
%   'AlphaThresholdValue' - if percentile: 1-100 (default: 50); 
%                           if absolute: threshold in same units as data (e.g., 5 for 5 mcV)
%   'SMIThresholdType' - 'auto' (default) or 'absolute'
%   'SMIThresholdValue' - if 'absolute': threshold in same units as data (e.g., 5 for 5 mcV)
%                        (not used if 'auto' is selected)
%
% Outputs:
%   AMI - Alpha Modulation Index
%   SMI - Slow Modulation Index

% Parse optional inputs
p = inputParser;
addParameter(p, 'AlphaThresholdType', 'percentile', @ischar);
addParameter(p, 'AlphaThresholdValue', 50, @isnumeric);
addParameter(p, 'SMIThresholdType', 'auto', @ischar);
addParameter(p, 'SMIThresholdValue', 5, @isnumeric);
parse(p, varargin{:});

alphaThresholdType = p.Results.AlphaThresholdType;
alphaThresholdValue = p.Results.AlphaThresholdValue;
smiThresholdType = p.Results.SMIThresholdType;
smiThresholdValue = p.Results.SMIThresholdValue;

% Define frequency bands
alpha_L = 8;
alpha_H = 14;
slow_L = 0.3;
slow_H = 4;

% Check input data format
if ~isvector(data)
    error('Input data must be a vector.');
end

% Ensure data is a column vector
if size(data, 1) == 1
    data = data';
end

% Filter signals to extract alpha and slow components
alpha = filter_butterworth(data, fs, alpha_L, alpha_H, 2);
vslow = filter_butterworth(data, fs, slow_L, slow_H, 2);

% Set processing boundaries
AMI_START = 1;  % MATLAB indexing starts at 1 (not 0)
AMI_END = length(data);
xlim_beg = 1;   % MATLAB indexing
xlim_end = length(data);

% Kalman filter parameters
cutoff_period = 5*60;
K = 2*pi/(fs * cutoff_period);
r_kalman = 1;
q_kalman = r_kalman * K^2 / (1-K); % variance def
r_kalman_slow = 10*r_kalman; % variance def

%% Compute AMI
beg = floor(threshold_TIME_L*fs) + 1;  % MATLAB indexing
endIdx = ceil(threshold_TIME_H*fs);

% Boundary check
if beg < 1 || endIdx > length(data)
    error('Threshold times out of data range. Adjust threshold_TIME_L and threshold_TIME_H.');
end

temp_signal = alpha(beg:endIdx);

% Determine threshold based on method
if strcmpi(alphaThresholdType, 'percentile')
    % Find alpha peaks in calibration window
    alpha_peaks = [];
    for i = 2:length(temp_signal)-1
        if (temp_signal(i) > temp_signal(i-1) && temp_signal(i) > temp_signal(i+1))
            alpha_peaks = [alpha_peaks, temp_signal(i)];
        end
    end
    
    % Check if there are any peaks
    if isempty(alpha_peaks)
        error('No alpha peaks detected in calibration window. Try a different window or consider using an absolute threshold.');
    end
    
    threshold = prctile(alpha_peaks, alphaThresholdValue);
    disp(['Using percentile threshold for AMI: ', num2str(threshold)]);
elseif strcmpi(alphaThresholdType, 'absolute')
    threshold = alphaThresholdValue;
    disp(['Using absolute threshold for AMI: ', num2str(threshold)]);
else
    error('Invalid threshold type. Use "percentile" or "absolute".');
end

% Threshold alpha signal
alpha_thresholded = (abs(alpha) >= threshold);
alpha_thresholded_padded = [0; alpha_thresholded];
alpha_crossings = diff(alpha_thresholded_padded);

% Apply window to binary signal
Window_size = 18;
N = length(alpha);
alpha_binary = zeros(N, 1);
for i = 1:N
    beg = max([1, i-Window_size]);
    endIdx = min([N, i]);
    alpha_binary(i) = max(alpha_thresholded(beg:endIdx));
end

% Detect up and down states
UP_START = [];
UP = [];
DOWN_START = [1];  % MATLAB indexing
DOWN = [];
prev = 0;
SUM = 0;
PREV_DOWN = 0;
PREV_UP = 0;
last_UP_start = 1;  % MATLAB indexing
last_DOWN_start = N;

UP_signal = zeros(N, 1);
UP_noZero_signal = zeros(N, 1);
DOWN_signal = zeros(N, 1);
DOWN_noZero_signal = zeros(N, 1);

for i = 1:N
    current = alpha_binary(i);
    if (prev == current)
        SUM = SUM + 1/fs;
    else
        if (prev == 0)
            UP_START = [UP_START, i];
            DOWN = [DOWN, SUM];
            if (last_UP_start ~= 1)  % MATLAB indexing
                DOWN_noZero_signal(last_UP_start:last_DOWN_start) = DOWN_noZero_signal(last_UP_start:last_DOWN_start) + PREV_DOWN;
                for j = last_DOWN_start:i-1
                    DOWN_noZero_signal(j) = max([DOWN_noZero_signal(j), PREV_DOWN]);
                end
            end
            last_UP_start = i;
            DOWN_signal(last_DOWN_start:i-1) = DOWN_signal(last_DOWN_start:i-1) + SUM;
            PREV_DOWN = SUM;
            SUM = 0;
        end
        if (prev == 1)
            DOWN_START = [DOWN_START, i];
            UP = [UP, SUM];
            UP_noZero_signal(last_DOWN_start:last_UP_start) = UP_noZero_signal(last_DOWN_start:last_UP_start) + PREV_UP;
            for j = last_UP_start:i-1
                UP_noZero_signal(j) = max([UP_noZero_signal(j), PREV_UP]);
            end
            last_DOWN_start = i;
            UP_signal(last_UP_start:i-1) = UP_signal(last_UP_start:i-1) + SUM;
            PREV_UP = SUM;
            SUM = 0;
        end
    end
    prev = current;
end

% Process UP/DOWN values using Kalman filter
MARKER_UP = kalmanFilter(UP_noZero_signal, q_kalman, r_kalman, AMI_START);
MARKER_DOWN = kalmanFilter(DOWN_noZero_signal, q_kalman, r_kalman, AMI_START);
gamma = 0.5;
AMI = MARKER_UP.^gamma ./ (MARKER_DOWN.^gamma + MARKER_UP.^gamma);

%% Compute SMI
beg = floor((threshold_TIME_L+xlim_beg-1)*fs) + 1;  % MATLAB indexing
endIdx = ceil((threshold_TIME_H+xlim_beg-1)*fs);
    
temp_data = vslow(beg:endIdx);

% Determine quiescence threshold based on method
if strcmpi(smiThresholdType, 'auto')
    max_level = 0;
    max_threshold = 0;
    RES = 1;
    beg_range = floor(min(temp_data)*0.9)*RES;
    end_range = floor(max(temp_data)*0.9)*RES;
    
    % Ensure we have a valid range
    if beg_range >= end_range
        beg_range = -10;
        end_range = 10;
    end
    
    all_levels = [];
    all_thresholds = [];

    % Find optimal threshold for detecting slow oscillations
    for i = beg_range:end_range
        temp_threshold = (1.0*i)/RES;
        level = computeRate(temp_data, temp_threshold, fs);
        all_levels = [all_levels, level];
        all_thresholds = [all_thresholds, temp_threshold];
        if (level > max_level)
            max_level = level;
            max_threshold = temp_threshold;
        end
    end

    % Check if we found any valid threshold
    if max_level == 0
        error('Could not find optimal slow oscillation threshold. Try using an absolute threshold.');
    end

    twoThird_level = 0.8*max_level;
    twoThird_threshold = max_threshold;
    
    % Find threshold where rate drops to 80% of max
    for i = 1:length(all_thresholds)
        if (all_thresholds(i) > max_threshold)
            if (all_levels(i) <= twoThird_level)
                twoThird_threshold = all_thresholds(i);
                break;
            end
        end
    end

    threshold = max_threshold;
    THE_THRESHOLD = threshold;
    quiescence_threshold = 0.5*twoThird_threshold + 0.5*THE_THRESHOLD;
    
    disp(['Using auto-computed threshold for SMI: ', num2str(quiescence_threshold)]);
elseif strcmpi(smiThresholdType, 'absolute')
    quiescence_threshold = smiThresholdValue;
    threshold = 0; % A value below the quiescence threshold to detect crossings
    THE_THRESHOLD = threshold;
    
    disp(['Using absolute threshold for SMI: ', num2str(quiescence_threshold)]);
else
    error('Invalid SMI threshold type. Use "auto" or "absolute".');
end
    
crossing_threshold = threshold;

% Process slow wave signal
data_thresholded = (vslow > crossing_threshold);
data_diff = diff([0; data_thresholded]);
rise = find(data_diff == 1);
fall = find(data_diff == -1);

if (~isempty(fall) && ~isempty(rise) && fall(1) < rise(1))
    fall = fall(2:end);
end

% Safety check for empty arrays
if isempty(rise) || isempty(fall)
    warning('No slow oscillation cycles detected. SMI will be unreliable.');
    SMI = zeros(size(data));
    return;
end

Nb_Segments = length(fall);

min_cycle_length = 1/slow_H;
vslow_crossings = zeros(size(vslow));
vslow_quiescence = zeros(size(vslow));

for s = 1:Nb_Segments-1
    if s <= length(rise) && rise(s) <= length(vslow) && s <= length(fall) && fall(s) <= length(vslow)
        start = rise(s);
        mid = fall(s);
        if s+1 <= length(rise)
            endIdx = rise(s+1);
            if start < mid && mid < endIdx && endIdx <= length(vslow)
                up = vslow(start:mid);
                [amplitude, peak_index] = max(up);
                vslow_crossings(start) = amplitude;
                if (amplitude <= quiescence_threshold)
                    vslow_quiescence(start:endIdx) = 1;
                end
            end
        end
    end
end
   
UP_START_vslow = find(vslow_crossings > 0);
SLOW_SIGNAL = zeros(size(vslow)) + 0.01;
SLOW_SIGNAL_DURATION = zeros(size(vslow)) + 0.01;
prev = 0;

for i = 1:length(UP_START_vslow)-1
    beg = UP_START_vslow(i);
    endIdx = UP_START_vslow(i+1);
    SLOW_SIGNAL(beg:endIdx) = SLOW_SIGNAL(beg:endIdx) + fs/(endIdx-beg);
    for j = beg:endIdx
        SLOW_SIGNAL_DURATION(j) = max([prev, (j-beg)/fs]);
    end
    prev = (endIdx-beg)/fs;
end

% Process quiescence periods
vslow_quiescence = (vslow <= quiescence_threshold);
diff_vslow_quiescence = diff([0; vslow_quiescence]);
UP_QUIESCENCE_vslow = find(diff_vslow_quiescence > 0);
DOWN_QUIESCENCE_vslow = find(diff_vslow_quiescence < 0);

if (~isempty(UP_QUIESCENCE_vslow) && ~isempty(DOWN_QUIESCENCE_vslow) && ...
        UP_QUIESCENCE_vslow(1) > DOWN_QUIESCENCE_vslow(1))
    DOWN_QUIESCENCE_vslow = DOWN_QUIESCENCE_vslow(2:end);
end

QUIESCENCE_SIGNAL = zeros(size(vslow));
QUIESCENCE_SIGNAL_REALTIME = zeros(size(vslow));
ONLY_QUIESCENCE_SIGNAL = zeros(size(vslow));
prev = 0;

for i = 1:length(UP_QUIESCENCE_vslow)-1
    if i <= length(DOWN_QUIESCENCE_vslow)
        beg = UP_QUIESCENCE_vslow(i);
        endIdx = DOWN_QUIESCENCE_vslow(i);
        if i+1 <= length(UP_QUIESCENCE_vslow)
            nex_beg = UP_QUIESCENCE_vslow(i+1);
            val = (endIdx-beg)/fs;
            
            for j = beg:endIdx
                QUIESCENCE_SIGNAL_REALTIME(j) = max([prev, (j-beg)/fs]);
            end
            
            if endIdx < nex_beg
                QUIESCENCE_SIGNAL_REALTIME(endIdx:nex_beg) = val;    
                QUIESCENCE_SIGNAL(beg:nex_beg) = val;
                ONLY_QUIESCENCE_SIGNAL(beg:endIdx) = val;
            end
            prev = val;
        end
    end
end
    
% Calculate SMI using Kalman filtered signals
% MARKER_SLOW = kalmanFilter(1./SLOW_SIGNAL_DURATION, q_kalman, r_kalman_slow, AMI_START);
inv_duration = 1 ./ max(SLOW_SIGNAL_DURATION, 1e-6);
MARKER_SLOW = kalmanFilter(inv_duration, q_kalman, r_kalman_slow, AMI_START);
MARKER_QUIESCENCE = kalmanFilter(QUIESCENCE_SIGNAL_REALTIME, q_kalman, r_kalman_slow, AMI_START);

beg = floor((threshold_TIME_L+xlim_beg-1)*fs) + 1;
endIdx = ceil((threshold_TIME_H+xlim_beg-1)*fs);
QUIESCENCE_LEVEL = mean(QUIESCENCE_SIGNAL(beg:endIdx));
MARKER_LEVEL = mean(MARKER_SLOW(beg:endIdx));

% Avoid division by zero
if QUIESCENCE_LEVEL == 0
    QUIESCENCE_LEVEL = eps;
end
if MARKER_LEVEL == 0
    MARKER_LEVEL = eps;
end

gamma = 2;
NUM = (MARKER_SLOW./MARKER_LEVEL).^gamma;
DENOM = ((MARKER_QUIESCENCE./QUIESCENCE_LEVEL)).^gamma;
SMI = NUM./(NUM+DENOM);

end

%% Helper Functions

function filtered = filter_butterworth(y, fs, lowcut, highcut, order)
    % Butterworth bandpass filter
    nyq = 0.5 * fs;
    low = lowcut / nyq;
    high = highcut / nyq;
    
    if low == 0
        [b, a] = butter(order, high, 'low');
    elseif isinf(high)
        [b, a] = butter(order, low, 'high');
    else
        [b, a] = butter(order, [low, high]);
    end
    
    % Filter in sections to avoid the filtfilt compatibility issue
    n = length(y);
    filtered = zeros(size(y));
    
    % Forward pass
    forward = filter(b, a, y);
    
    % Backward pass (to get zero-phase)
    backward = flip(filter(b, a, flip(forward)));
    
    filtered = backward;
end

function [x, p] = kalmanInitCondition(z, q, r, START)
% Initialize Kalman filter with reversed signal
z_rev = flipud(z);
[x_rev, p_rev] = kalmanFilter_core(z_rev, q, r, 1, 1, 1);
x = flipud(x_rev);
p = flipud(p_rev);
x = x(START);
p = p(START);
end

function [x, p] = kalmanFilter_core(z, q, r, START, x0, p0)
% Core Kalman filter implementation
N = length(z);
x = zeros(N, 1);
p = zeros(N, 1);
x(1:START) = x0;
p(1:START) = p0;

for i = START+1:N
    % predict
    x(i) = x(i-1);
    p(i) = p(i-1) + q;
    % update
    K = p(i)/(p(i) + r);
    x(i) = (1-K)*x(i) + K*z(i);
    p(i) = (1-K)*p(i);
end
end

function [xs, ps, x, p, k, l] = kalmanSmoother_core(z, q, r, START, x0, p0, END)
% Core Kalman smoother implementation
N = length(z);
x = zeros(N, 1);
p = zeros(N, 1);
k = zeros(N, 1);
x(1:START) = x0;
p(1:START) = p0;

for i = START+1:N
    % predict
    x(i) = x(i-1);
    p(i) = p(i-1) + q;
    % update
    K = p(i)/(p(i) + r);
    x(i) = (1-K)*x(i) + K*z(i);
    p(i) = (1-K)*p(i);
    k(i) = K;
end

xs = zeros(N, 1);
ps = zeros(N, 1);
l = zeros(N, 1);
xs(END) = x(END);
ps(END) = p(END);

for i = START:END-1
    j = END + START - i;
    L = p(j)/p(j+1);
    xs(j) = x(j) + L*(xs(j+1) - x(j+1));
    ps(j) = p(j) + L*L*(ps(j+1) - p(j+1));
    l(j) = L;
end
end

function [x, p] = kalmanFilter(z, q, r, START)
% Kalman filter implementation
[x0, p0] = kalmanInitCondition(z, q, r, START);
[x, p] = kalmanFilter_core(z, q, r, START, x0, p0);
end

function [xs, ps, x, p, k, l] = kalmanSmoother(z, q, r, START, END)
% Kalman smoother implementation
[x0, p0] = kalmanInitCondition(z, q, r, START);
[xs, ps, x, p, k, l] = kalmanSmoother_core(z, q, r, START, x0, p0, END);
end

function rate = computeRate(vslow_data, crossing_threshold, fs)
% Compute rate of slow oscillation cycles
data_thresholded = (vslow_data > crossing_threshold);
data_diff = diff([0; data_thresholded]);
rise = find(data_diff == 1);
fall = find(data_diff == -1);

if ~isempty(fall) && ~isempty(rise) && (fall(1) < rise(1))
    fall = fall(2:end);
end

Nb_Segments = length(fall);
if Nb_Segments == 0
    rate = 0;
    return;
end

vslow_crossings = zeros(size(vslow_data));
for s = 1:Nb_Segments-1
    if s <= length(rise) && s+1 <= length(rise)
        start = rise(s);
        mid = fall(s);
        endIdx = rise(s+1);
        if start <= length(vslow_data) && mid <= length(vslow_data) && endIdx <= length(vslow_data)
            up = vslow_data(start:mid);
            [amplitude, peak_index] = max(up);
            vslow_crossings(start) = amplitude;
        end
    end
end

UP_START_vslow = find(vslow_crossings > 0);
SLOW_SIGNAL = zeros(size(vslow_data));

for i = 1:length(UP_START_vslow)-1
    beg = UP_START_vslow(i);
    endIdx = UP_START_vslow(i+1);
    if beg < endIdx && endIdx <= length(SLOW_SIGNAL)
        SLOW_SIGNAL(beg:endIdx-1) = fs/(endIdx-beg);
    end
end

rate = mean(SLOW_SIGNAL);
end