%% Level-Crossing vs Uniform-Sample Autozero: thermal + flicker, Welch PSD
clear; clc;
close all;

%% ------------------------ User parameters -----------------------------
Vamp    = 6;               % sine amplitude [V]
fin     = 1;               % sine frequency [Hz]
LSB     = 500e-3;          % quantization step [V] (spike step)

Fs      = 1e6;             % simulation sample rate [Hz]
Tsec    = 10;              % total duration [s]
rng(42)                    % reproducibility

% Welch PSD controls (CAUTION: huge windows at Fs=1e6)
windowTime  = 10;          % seconds per Welch segment
windowSize  = round(Fs * windowTime);
overlap     = round(windowSize/2);

% Full-scale definition (for dBFS)
FS_Vpp   = 6;
FS_RMS   = FS_Vpp / (2 * sqrt(2));
FS_Power = FS_RMS^2;
FS_dB    = 10*log10(FS_Power); %#ok<NASGU>

% Noise-band & guards
flo       = 1;             % Hz, lower bound for noise integration
fbw       = 500;           % Hz, upper bound for noise integration
harmonics = 10;            % number of harmonics to exclude from noise
%% ----------------------------------------------------------------------

% ---- Thermal noise (white, input-referred) ----
addThermalNoise = true;
noiseRMS        = 0.1;     % [V_RMS] (set to 0 to disable)

% ---- Flicker noise (1/f, input-referred) ----
addFlickerNoise = true;
flickerASD_1Hz  = 50e-3;    % [V/sqrt(Hz)] @ 1 Hz (magnitude knob)

% ---- Auto-zero settings ----
useUniformAZ     = true;   % uniform-sampling AZ (NO LC)
Fsample_Hz       = 1e3;    % uniform sampling + AZ rate [Hz]
useEventAZ       = true;   % event-triggered AZ on LC path (YES LC)
%% ----------------------------------------------------------------------

%% ------------------------ Timebase and stimulus ------------------------
N   = round(Tsec*Fs);
t   = (0:N-1).'/Fs;
x0  = Vamp * sin(2*pi*fin*t);            % clean input (no noise yet)

% Build noises separately so AZ only cancels flicker (ideal behavior)
nth = zeros(N,1);
if addThermalNoise && noiseRMS>0
    nth = noiseRMS * randn(N,1);         % thermal (white)
end

nfk = zeros(N,1);
if addFlickerNoise && flickerASD_1Hz>0
    nfk = make_flicker_noise(N, Fs, Tsec, flickerASD_1Hz);  % memory-light
end

%% ---------------- Paths & reconstructions ------------------------------
% A) Level-crossing, NO autozero (baseline)
[spk_on0, spk_off0, y_rec_noAZ] = lc_encode_eventAZ(x0, LSB, nth, nfk, 'none');

% B) Level-crossing, EVENT-based ideal autozero (flicker only)
if useEventAZ
    [spk_onE, spk_offE, y_rec_eAZ] = lc_encode_eventAZ(x0, LSB, nth, nfk, 'event');
end

% C) UNIFORM sampling (+AZ at the same rate) with ZOH back to Fs grid (NO LC)
if useUniformAZ
    y_rec_uAZ = uniform_sample_autozero_ZOH(x0, nth, nfk, Fs, Fsample_Hz);
end

%% ---------------- PSD + metrics (Welch) --------------------------------
[f0, PSD0_dBFSbin, M0] = your_welch_psd_and_metrics(y_rec_noAZ, Fs, fin, ...
    windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime);

if useEventAZ
    [fE, PSDE_dBFSbin, ME] = your_welch_psd_and_metrics(y_rec_eAZ, Fs, fin, ...
        windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime);
end
if useUniformAZ
    [fU, PSDU_dBFSbin, MU] = your_welch_psd_and_metrics(y_rec_uAZ, Fs, fin, ...
        windowSize, overlap, FS_Power, flo, fbw, harmonics, windowTime);
end

fprintf('\n=== Metrics (LC: No AZ) ===\n');
fprintf('SNR   : %6.2f dB\n',  M0.SNR_dB);
fprintf('SNDR  : %6.2f dB\n',  M0.SNDR_dB);
fprintf('SFDR  : %6.2f dB\n',  M0.SFDR_dB);
fprintf('ENOB  : %6.3f bits\n',M0.ENOB_bits);

if useEventAZ
    fprintf('\n=== Metrics (LC: Event-based AZ) ===\n');
    fprintf('SNR   : %6.2f dB\n',  ME.SNR_dB);
    fprintf('SNDR  : %6.2f dB\n',  ME.SNDR_dB);
    fprintf('SFDR  : %6.2f dB\n',  ME.SFDR_dB);
    fprintf('ENOB  : %6.3f bits\n',ME.ENOB_bits);
end

if useUniformAZ
    fprintf('\n=== Metrics (Uniform sample + AZ @ %g Hz) ===\n', Fsample_Hz);
    fprintf('SNR   : %6.2f dB\n',  MU.SNR_dB);
    fprintf('SNDR  : %6.2f dB\n',  MU.SNDR_dB);
    fprintf('SFDR  : %6.2f dB\n',  MU.SFDR_dB);
    fprintf('ENOB  : %6.3f bits\n',MU.ENOB_bits);
end

% Spike rate sanity (events/second) for LC paths
rate_noAZ = (sum(double(spk_on0))+sum(double(spk_off0)))/Tsec;
fprintf('\nSpike rate (LC, no AZ): %.2f events/s\n', rate_noAZ);
if useEventAZ
    rate_eAZ  = (sum(double(spk_onE))+sum(double(spk_offE)))/Tsec;
    fprintf('Spike rate (LC, event AZ): %.2f events/s\n', rate_eAZ);
end
if useUniformAZ
    fprintf('Uniform path: no spikes (sampled @ %.2f Hz)\n', Fsample_Hz);
end

%% ------------------------------ Plots ----------------------------------
% Time-domain: input + three reconstructions
tview = t < Tsec;
figure('Name','Time domain');
subplot(3,1,1);
plot(t(tview), x0(tview) + nth(tview) + nfk(tview), 'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('x(t) [V]');
title('Input sine + thermal + flicker');

subplot(3,1,2);
plot(t(tview), y_rec_noAZ(tview), 'LineWidth',1); grid on;
xlabel('Time [s]'); ylabel('\itŷ\rm [V]'); title('LC Reconstruction (No AZ)');

subplot(3,1,3);
hold on; grid on;
if useUniformAZ, plot(t(tview), y_rec_uAZ(tview), 'LineWidth',1); end
if useEventAZ,   plot(t(tview), y_rec_eAZ(tview), 'LineWidth',1); end
xlabel('Time [s]'); ylabel('\itŷ\rm [V]');
lg = {};
if useUniformAZ, lg{end+1} = sprintf('Uniform sample + AZ @ %g Hz',Fsample_Hz); end %#ok<AGROW>
if useEventAZ,   lg{end+1} = 'LC: Event-based AZ'; end %#ok<AGROW>
legend(lg,'Location','best'); title('Reconstructions with ideal autozero');

% Spike rasters (zoomed; LC baseline for readability)
figure('Name','Spike trains (zoomed, LC No AZ)');
Tzoom = min(1e-3, Tsec);
tzoom = t < Tzoom;
yyaxis left
stem(t(tzoom)*1e3, spk_on0(tzoom), '.', 'Marker','none'); hold on;
ylabel('ON spikes')
yyaxis right
stem(t(tzoom)*1e3, -spk_off0(tzoom), '.', 'Marker','none');
ylabel('OFF spikes (negative)');
grid on; xlabel('Time [ms]');
title('ON/OFF spike trains (LC, No AZ)');

% Welch PSD plots (log frequency)
figure('Name','Welch PSD (log f) - AZ comparison');
k0 = find(f0>0,1,'first');
semilogx(f0(k0:end), PSD0_dBFSbin(k0:end), 'LineWidth',1); hold on; grid on; grid minor;
leg = {'LC: No AZ'};
if useEventAZ
    kE = find(fE>0,1,'first');
    semilogx(fE(kE:end), PSDE_dBFSbin(kE:end), 'LineWidth',1);
    leg{end+1} = 'LC: Event-based AZ';
end
if useUniformAZ
    kU = find(fU>0,1,'first');
    semilogx(fU(kU:end), PSDU_dBFSbin(kU:end), 'LineWidth',1);
    leg{end+1} = sprintf('Uniform sample + AZ @ %g Hz',Fsample_Hz);
end
xlabel('Frequency [Hz]');
ylabel('PSD [dBFS]');
legend(leg,'Location','southwest');
%title(sprintf('Welch PSD (window %.3f s, overlap %d samples, N=%.0f)', ...
%    windowTime, overlap, N));

if exist('micasplot','file'), micasplot, end
width = 800; height = 400;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['./az_model.png']);

%% -------- UPDATED: Bar plot comparing SNDR, SFDR, SNR (no ENOB) --------
metricCats = categorical({'SNR','SNDR','SFDR'});
metricCats = reordercats(metricCats, {'SNR','SNDR','SFDR'});

Y = [];                 % 3 x numCases matrix (rows: metrics; cols: cases)
caseLabels = {};

% LC: No AZ (always present)
Y = [Y, [M0.SNR_dB; M0.SNDR_dB; M0.SFDR_dB]];
caseLabels{end+1} = 'LC: No AZ';

% LC: Event-based AZ (optional)
if useEventAZ
    Y = [Y, [ME.SNR_dB; ME.SNDR_dB; ME.SFDR_dB]];
    caseLabels{end+1} = 'LC: Event AZ';
end

% Uniform sample + AZ (optional)
if useUniformAZ
    Y = [Y, [MU.SNR_dB; MU.SNDR_dB; MU.SFDR_dB]];
    caseLabels{end+1} = sprintf('Uniform + AZ @ %g Hz', Fsample_Hz);
end

figure('Name','Metric comparison (bar) – SNR/SNDR/SFDR');
bar(metricCats, Y, 'grouped'); grid on;
ylabel('dB');
legend(caseLabels,'Location','best');
% title('SNR / SNDR / SFDR across methods');

if exist('micasplot','file'), micasplot, end
width = 800; height = 300;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, './az_metrics_bar.png');

%% =========================== FUNCTIONS ================================
function [spk_on, spk_off, yhat] = lc_encode_eventAZ(x0, LSB, nth, nfk, mode)
% Level-crossing encoder with optional ideal EVENT-based autozero of FLICKER.
% mode = 'none' | 'event'
    N = numel(x0);
    spk_on  = zeros(N,1,'int16');
    spk_off = zeros(N,1,'int16');
    yhat    = zeros(N,1);
    ycurr   = 0;
    half    = LSB/2;

    az_state = 0;   % last sampled flicker offset (for event-AZ)

    for n = 1:N
        % Comparator input with ideal flicker cancellation (thermal unaffected)
        if strcmp(mode,'none')
            xi = x0(n) + nth(n) + nfk(n);
        else
            xi = x0(n) + nth(n) + (nfk(n) - az_state);
        end

        e = xi - ycurr;
        crossed = false;

        while e >= half
            spk_on(n) = spk_on(n) + 1;
            ycurr = ycurr + LSB;
            e = xi - ycurr;
            crossed = true;
        end
        while e <= -half
            spk_off(n) = spk_off(n) + 1;
            ycurr = ycurr - LSB;
            e = xi - ycurr;
            crossed = true;
        end

        % Ideal event-based AZ: update after any crossing
        if crossed && strcmp(mode,'event')
            az_state = nfk(n);
        end

        yhat(n) = ycurr;
    end
end

function y = uniform_sample_autozero_ZOH(x0, nth, nfk, Fs, Fsample_Hz)
% Uniformly sample at Fsample_Hz; at each sample, ideal-AZ cancels flicker only.
% Then zero-order hold (ZOH) to the high-rate grid (Fs) so PSDs are comparable.
    N  = numel(x0);
    Ns = max(1, round(Fs / Fsample_Hz));     % samples per hold at Fs
    idx = 1:Ns:N;                             % sample indices on Fs grid
    y   = zeros(N,1);

    for k = 1:numel(idx)
        n    = idx(k);
        az   = nfk(n);                        % ideal measurement of flicker
        x_s  = x0(n) + nth(n) + (nfk(n) - az); % thermal left, flicker canceled
        n1   = n;
        n2   = min(N, n + Ns - 1);
        y(n1:n2) = x_s;                       % hold until next sample
    end
end

function nfk = make_flicker_noise(N, Fs, Tsec, asd1Hz)
% Memory-light 1/f noise via sum of first-order lowpasses (octave sum).
% Scales output so total variance ~ asd^2 * ln(fmax/fmin), giving PSD~1/f.
    if asd1Hz<=0
        nfk = zeros(N,1);
        return;
    end
    fmin = max(1/Tsec, 1/(N/Fs));   % avoid too-small low corner
    fmax = Fs/2;

    % Choose log-spaced poles across the band (tune K for speed vs. fidelity)
    K  = 16;
    fp = logspace(log10(fmin), log10(fmax), K).';

    y = zeros(N,1);
    for k = 1:K
        a = exp(-2*pi*fp(k)/Fs);           % AR(1) pole
        b = sqrt(max(1 - a^2, 0));         % stable gain term
        w = randn(N,1);
        z = filter(b, [1 -a], w);          % first-order colored process
        y = y + z ./ sqrt(fp(k));          % ~1/sqrt(f) weighting => PSD ~1/f
    end

    % DC removal & scale to target variance ~ asd^2 * ln(fmax/fmin)
    y  = y - mean(y);
    sx = std(y);
    targetVar = (asd1Hz^2) * max(log(fmax/fmin), 0);
    if sx>0 && targetVar>0
        y = y * (sqrt(targetVar)/sx);
    else
        y = zeros(size(y));
    end
    nfk = y;
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

