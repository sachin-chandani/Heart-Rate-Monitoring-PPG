load './DATA_06_TYPE02.mat';
load 'DATA_06_TYPE02_BPMtrace.mat';
data = sig;

%% Initialization
I = data(2:end, :)';
I_sz = size(I);
Fs = 125;
fmin = 1.2;
fmax = 2.5;
tmax = floor(I_sz(1) / Fs);
flag = 0;
if rem(I_sz(1), Fs) == 0
    flag = 1;
end
delt = 8;
step = 2;
s = 0:2:(tmax - delt);
s = (s * Fs) + 1;
s(end) = s(end) - flag;
dist = delt * Fs;
dilation = 64 * Fs * 16;

%% BandPass Filter
I(:, 1) = bandpass(I(:, 1), [0.8 3.0], Fs);
I(:, 2) = bandpass(I(:, 2), [0.8 3.0], Fs);

%% Calculate acceleration
accx = I(:, 3) .^ 2;
accy = I(:, 4) .^ 2;
accz = I(:, 5) .^ 2;
acc = accx + accy + accz;

%% Denoising
I(:, 1) = filter_noise(I(:, 1), accx);
I(:, 2) = filter_noise(I(:, 2), accx);

I(:, 1) = filter_noise(I(:, 1), accy);
I(:, 2) = filter_noise(I(:, 2), accy);

I(:, 1) = filter_noise(I(:, 1), accz);
I(:, 2) = filter_noise(I(:, 2), accz);

% Denoising step by step
filter_delt = 30;
filter_step = 10;
st = 0:filter_step:(tmax - filter_delt);
st = st * Fs + 1;
filter_dist = filter_delt * Fs;
for i = st
    ind = i:i+filter_dist;
    I(ind, 1) = filter_noise(I(ind, 1), accx(ind));
    I(ind, 2) = filter_noise(I(ind, 2), accx(ind));
    
    I(ind, 1) = filter_noise(I(ind, 1), accy(ind));
    I(ind, 2) = filter_noise(I(ind, 2), accy(ind));
    
    I(ind, 1) = filter_noise(I(ind, 1), accz(ind));
    I(ind, 2) = filter_noise(I(ind, 2), accz(ind));

end

%% HR Estimation
Y1 = [];
Y2 = [];
win = 7;
win1 = win;
win2 = win;

for i = s
    PPG1 = I(i:(i + dist), 1);
    PPG2 = I(i:(i + dist), 2);
    PPG1 = padarray(PPG1, dilation - dist, 'post');
    PPG2 = padarray(PPG2, dilation - dist, 'post');
    FFTPPG1 = fft(PPG1);
    FFTPPG2 = fft(PPG2);
    
    first = ceil(fmin * dilation / Fs);
    last = ceil(fmax * dilation / Fs);
    
    temp = FFTPPG1(first:last, :);
    temp = abs(temp);
    [~, peak_ppg1] = max(temp);
    peak_ppg1 = peak_ppg1 + first - 1;
    
    temp = FFTPPG2(first:last, :);
    temp = abs(temp);
    [~, peak_ppg2] = max(temp);
    peak_ppg2 = peak_ppg2 + first - 1;

    fppg1 = peak_ppg1 * Fs / dilation;
    fppg2 = peak_ppg2 * Fs / dilation;

    HR1 = fppg1 * 60;
    HR2 = fppg2 * 60;
    
    a = 0.999;
    disc = 0.8;
    if i ~= 1
        if abs(HR1 - Y1(end)) > win1
            HR1 = (1 - a)*HR1 + a*Y1(end);
        else
            win1 = win;
        end
        if abs(HR2 - Y2(end)) > win2
            HR2 = (1 - a)*HR2 + a*Y2(end);
        else
            win2 = win;
        end
        win1 = win1 * disc;
        win2 = win2 * disc;
    end
    Y1 = [Y1; HR1];
    Y2 = [Y2; HR2];
end

%% HR Smoothening
L = 9;
Yavg = max(Y1, Y2);

%% Error calculation
err_hr1 = Y1 - BPM0;
err_hr1 = err_hr1 .^ 2;
sum(err_hr1)

err_hr2 = Y2 - BPM0;
err_hr2 = err_hr2 .^ 2;
sum(err_hr2)

err_hravg = Yavg - BPM0;
err_hravg = abs(err_hravg);
sz = size(err_hravg);
sum(err_hravg) / sz(1)

%% Functions
% Noise filter
function [Y] = filter_noise(I, N)
    M      = 7;
    delta  = 0.1;
    P0     = (1/delta)*eye(M,M);
    rlsfilt = dsp.RLSFilter(M,'InitialInverseCovariance',P0);
    [y e] = rlsfilt(N, I);
    Y = e;
end

function [ind] = find_peak(I, offset, base, Fs, dilation, delta)
    flag = true;
    while flag
        [~, loc] = max(I);
        ind = loc + offset - 1;
        hr = ind * Fs / dilation * 60;
        if abs(hr - base) <= delta
            flag = false;
        else
            I(loc) = 0;
        end
    end
end