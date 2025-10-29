%% Level-Crossing (ON/OFF spikes) with Welch PSD (your method) + random blanking
clear; clc;
% close all;

%% ------------------------ User parameters -----------------------------
Vamp    = 6;              % sine amplitude [V]
fin     = 1;              % sine frequency [Hz]  (a.k.a. freqin)
LSB     = 500e-3;         % quantization step [V]  (spike step)

Fs      = 1e6;            % simulation sample rate [Hz] (samplingRate)
Tsec    = 10;             % total duration [s]
rng(42)                   % reproducible blanking pattern

% Welch PSD controls
windowTime  = 10;          % seconds per Welch segment
windowSize  = round(Fs * windowTime);
overlap     = round(windowSize/2);

% Full-scale definition (for dBFS)
FS_Vpp  = 6;
FS_RMS  = FS_Vpp / (2 * sqrt(2));
FS_Power = FS_RMS^2;
FS_dB   = 10*log10(FS_Power); %#ok<NASGU>

% Noise-band & guards
flo = 3;                  % Hz, lower bound for noise integration
fbw = 500;                % Hz, upper bound for noise integration
harmonics = 10;           % number of harmonics to exclude from noise
%% ----------------------------------------------------------------------

% Blanking spec
blank_us    = 150;        % each blanking window length [microseconds]
gapMin_ms   = 0.1/fin;        % minimum inter-blank spacing [ms]
gapMax_ms   = 10/fin;         % maximum inter-blank spacing [ms]

% Optional thermal noise injection
addThermalNoise = true;   % <-- toggle this flag
noiseRMS = 0.1;          % noise RMS amplitude [V] (set to 0 for no noise)
%% ----------------------------------------------------------------------

%% ------------------------ Timebase and stimulus ------------------------
N   = round(Tsec*Fs);
t   = (0:N-1).'/Fs;
x   = Vamp * sin(2*pi*fin*t);

if addThermalNoise
    noise = noiseRMS * randn(size(x));
    x = x + noise;
end

%% ------------------------ Encode to ON/OFF spikes ----------------------
[spk_on, spk_off, y_rec_ideal] = lc_encode_and_reconstruct(x, LSB);

%% ---------------- PSD + metrics (ideal reconstruction, Welch) ----------
[f1, PSD1_dBFSbin, M1] = your_welch_psd_and_metrics(y_rec_ideal, Fs, fin, ...
    windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime);

fprintf('\n=== Metrics (no blanking) ===\n');
fprintf('SNR   : %6.2f dB\n',  M1.SNR_dB);
fprintf('SNDR  : %6.2f dB\n',  M1.SNDR_dB);
fprintf('SFDR  : %6.2f dB\n',  M1.SFDR_dB);
fprintf('ENOB  : %6.3f bits\n',M1.ENOB_bits);

%% -------------------- Random blanking of spike trains ------------------
mask_blank = make_blanking_mask(N, Fs, blank_us, gapMin_ms, gapMax_ms);
spk_on_b   = spk_on; spk_off_b = spk_off;
spk_on_b(mask_blank)  = 0;
spk_off_b(mask_blank) = 0;

% Reconstruction with induced error (explicit double to avoid integer ops)
y_rec_err = LSB * (cumsum(double(spk_on_b)) - cumsum(double(spk_off_b)));

%% ------------- PSD + metrics (with blanking error, Welch) -------------
[f2, PSD2_dBFSbin, M2] = your_welch_psd_and_metrics(y_rec_err, Fs, fin, ...
    windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime);

fprintf('\n=== Metrics (with blanking) ===\n');
fprintf('SNR   : %6.2f dB\n',  M2.SNR_dB);
fprintf('SNDR  : %6.2f dB\n',  M2.SNDR_dB);
fprintf('SFDR  : %6.2f dB\n', M2.SFDR_dB);
fprintf('ENOB  : %6.3f bits\n',M2.ENOB_bits);

%% ------------------------------ Plots ----------------------------------
% Time-domain view
tview = t < Tsec;
figure('Name','Time domain');
subplot(3,1,1);
plot(t(tview), x(tview), 'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('x(t) [V]'); title('Input sine (+ optional noise)');

subplot(3,1,2);
plot(t(tview), y_rec_ideal(tview), 'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('\itŷ\rm (ideal) [V]');
title('Reconstruction from spikes (no blanking)');

subplot(3,1,3);
plot(t(tview), y_rec_err(tview), 'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('\itŷ\rm (blanked) [V]');
title('Reconstruction from spikes (with blanking)');

% Spike rasters (zoomed for readability)
figure('Name','Spike trains (zoomed)');
Tzoom = min(1e-3, Tsec);
tzoom = t < Tzoom;
yyaxis left
stem(t(tzoom)*1e3, spk_on(tzoom), '.', 'Marker','none'); hold on;
ylabel('ON spikes')
yyaxis right
stem(t(tzoom)*1e3, -spk_off(tzoom), '.', 'Marker','none');
ylabel('OFF spikes (negative)');
grid on; xlabel('Time [ms]');
title('ON/OFF spike trains (counts per sample)');

% Welch PSD plots (log frequency)
figure('Name','Welch PSD (log f)');
kpos1 = find(f1>0,1,'first');
kpos2 = find(f2>0,1,'first');
semilogx(f1(kpos1:end), PSD1_dBFSbin(kpos1:end), 'LineWidth',1); hold on; grid on; grid minor;
semilogx(f2(kpos2:end), PSD2_dBFSbin(kpos2:end), 'LineWidth',1);
xlabel('Frequency [Hz]');
ylabel('PSD per bin [dBFS]');
legend('No blanking','With blanking','Location','southwest');
title(sprintf('Welch PSD (window %.3f s, overlap %d samples, N=%.0f)', ...
    windowTime, overlap, N));

% Blanking mask
figure('Name','Blanking mask');
stairs(t*1e3, mask_blank, 'LineWidth',1); grid on;
xlabel('Time [ms]'); ylabel('Blanked? (1=yes)');
title('Random 150 \\mus blanking windows (spike trains zeroed)');

%% =========================== FUNCTIONS ================================

function [spk_on, spk_off, yhat] = lc_encode_and_reconstruct(x, LSB)
    N = numel(x);
    spk_on  = zeros(N,1,'int16');
    spk_off = zeros(N,1,'int16');
    yhat    = zeros(N,1);
    ycurr   = 0;
    half    = LSB/2;
    for n = 1:N
        e = x(n) - ycurr;
        while e >= half
            spk_on(n) = spk_on(n) + 1;
            ycurr = ycurr + LSB;
            e = x(n) - ycurr;
        end
        while e <= -half
            spk_off(n) = spk_off(n) + 1;
            ycurr = ycurr - LSB;
            e = x(n) - ycurr;
        end
        yhat(n) = ycurr;
    end
end

function mask = make_blanking_mask(N, Fs, blank_us, gapMin_ms, gapMax_ms)
    mask = false(N,1);
    blankSamp = max(1, round(blank_us*1e-6 * Fs));
    startIdx  = 1;
    while startIdx <= N
        gapSamp = round( (gapMin_ms + (gapMax_ms-gapMin_ms)*rand) * 1e-3 * Fs );
        startIdx = startIdx + gapSamp;
        if startIdx > N, break; end
        stopIdx = min(N, startIdx + blankSamp - 1);
        mask(startIdx:stopIdx) = true;
        startIdx = stopIdx + 1;
    end
end

function [f, PSD_dBFSbin, M] = your_welch_psd_and_metrics( ...
        sig, Fs, freqin, windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime)

    y = double(sig(:));
    y = y - mean(y);

    win = hann(windowSize,'periodic');
    noverlap = overlap;
    nfft = windowSize;
    [Pxx, f] = pwelch(y, win, noverlap, nfft, Fs, 'psd');  % one-sided

    deltaF = Fs / windowSize;
    windowGain = sum(win.^2) / windowSize;
    psdEstimate = (Pxx * deltaF) / windowGain;
    psdEstimate = psdEstimate ./ FS_Power;

    if ~isempty(freqin) && freqin > 0
        [~, fundIndex] = min(abs(f - freqin));
    else
        [~, fundIndex] = max(psdEstimate(2:end)); fundIndex = fundIndex + 1;
    end

    sideband = 2 / windowTime;
    signalIndices = find(f >= (f(fundIndex) - sideband) & f <= (f(fundIndex) + sideband));
    signalPower   = sum(psdEstimate(signalIndices));

    harmonicIndices = zeros(harmonics,1);
    for h = 2:harmonics+1
        f_h = h * f(fundIndex);
        if f_h > f(end), break; end
        [~, harmonicIndices(h-1)] = min(abs(f - f_h));
    end
    harmonicIndices = harmonicIndices(harmonicIndices>0);
    harmonicPower   = sum(psdEstimate(harmonicIndices));

    lowerFreqLimit = flo + 1/windowTime;
    upperFreqLimit = fbw;
    noiseIndices = find(f > lowerFreqLimit & f <= upperFreqLimit);
    noiseIndices = setdiff(noiseIndices, harmonicIndices);
    noiseIndices = setdiff(noiseIndices, signalIndices);
    noisePower = sum(psdEstimate(noiseIndices));

    distortionPower = harmonicPower;
    SNDR_dB = 10*log10(signalPower / max(noisePower + distortionPower, eps));
    SNR_dB  = 10*log10(signalPower / max(noisePower, eps));
    if isempty(noiseIndices)
        SFDR_dB = NaN;
    else
        spuriousPower = max(psdEstimate(noiseIndices));
        SFDR_dB = 10*log10(signalPower / max(spuriousPower, eps));
    end
    ENOB_bits = (SNDR_dB - 1.76)/6.02;

    PSD_dBFSbin = 10*log10(psdEstimate + eps);

    M = struct('SNR_dB', SNR_dB, 'SNDR_dB', SNDR_dB, ...
               'SFDR_dB', SFDR_dB, 'ENOB_bits', ENOB_bits, ...
               'fundHz', f(fundIndex), 'sidebandHz', sideband, ...
               'noiseBandHz', [lowerFreqLimit, upperFreqLimit], ...
               'deltaF_Hz', deltaF, 'windowSize', windowSize, 'overlap', overlap);
end
