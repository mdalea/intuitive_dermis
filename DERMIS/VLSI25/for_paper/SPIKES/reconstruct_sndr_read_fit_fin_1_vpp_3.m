
%close all
clear all

name = 'sample_5_adc_ibias_5n_ibias_adapt_150n_islew_1u_vtrig_0p5_vpp_3_fin_1_vddadc_3_maxspikes_50000_ackdelay_1us-2024-11-22_12-11-41_sorted'

tlo=10+0.15+0.98; %0; %0.32; %0.94; %1.88; %3.85; 
thi=tlo+100;%113.9; %600; %11.9; %20.5; %11.88; %14.66; %23.85; %53.6;
mr_cnt = 1; %4;
tduration = thi-tlo; %10; %20; %seconds
t=[0:1:tduration*1e6]; % 400ms; 1us period spacing
sig = zeros(mr_cnt,length(t)); 
spline_sig = zeros(size(sig));	   

freqin=1; % fundamental frequency, known from siggen
flo=1;
fbw=50; %300; %1e3; % signal bandwidth
samplingRate = 1e6; %for %500e3; % Sampling rate (5 MHz in this example)
windowTime = 10; %0.5; %0.1; %1; %change to 1?
windowSize = samplingRate * windowTime;
overlap = windowSize / 2;

frame_data_on=zeros(2,2); frame_data_off=zeros(2,2);

%==== CT_MODE = '1' (spikes reversed) =====%
delta_hi=-1.33*10e-3; delta_lo=-1.33*10e-3; %2.75*delta_hi %18e-3
%delta_hi=-0.00123; delta_lo=-0.00084;
%==== CT_MODE = '0'  =====%
%delta_hi=18e-3; delta_lo=18e-3; %2.75*delta_hi %18e-3

Vpp = 3000e-3; %200e-3; % 200mV peak-to-peak amplitude

% Generate reference sinusoidal signal
t_ref = (0:length(t)-1) / samplingRate; % Time vector
ref_signal = (Vpp/2) * sin(2 * pi * freqin * t_ref); % 200mV peak-to-peak sinusoidal

%% READ spike data
filename = ['out\',name,'\1.bs2']

figure
TD = Read_Ndataset_longTs(filename);
%Encode_Ndataset([folderPath,'/',num2str(60001)],TD); % for testdata, for just for datasort_timestep to work
TD=dataset_sort_timestep_arduino(['out\',name],1);
%TD=dataset_sort_timestep_arduino(folderPath,1);      



% filter spikes within time range
TD.x = TD.x(TD.ts > tlo*1e6);
TD.y = TD.y(TD.ts > tlo*1e6);
TD.p = TD.p(TD.ts > tlo*1e6);
TD.ts = TD.ts(TD.ts > tlo*1e6);

TD.x = TD.x(TD.ts < thi*1e6);
TD.y = TD.y(TD.ts < thi*1e6);
TD.p = TD.p(TD.ts < thi*1e6);
TD.ts = TD.ts(TD.ts < thi*1e6);

delta_lo = delta_hi*(length(find(TD.p==1))/length(find(TD.p==2)))

subplot(3,1,1);
hold all;

i=2
length(TD.ts)
TD.ts = TD.ts - min(TD.ts); %normalize spike time start to 0
idx_within_time = find(TD.ts<tduration*1e6);

sparse_sig=[0];
%for spikeCount = 1:length(TD.ts)
for spikeCount = 1:max(idx_within_time)

    TD.ts(spikeCount)/1e6
    if TD.p(spikeCount)==1  % ON spike
        plot([TD.ts(spikeCount)/1e6 TD.ts(spikeCount)/1e6], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),2)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),2)+0.4], 'k');
        if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
            frame_data_on(TD.y(spikeCount)+1,TD.x(spikeCount)+1) = frame_data_on(TD.y(spikeCount)+1,TD.x(spikeCount)+1) + 1; % accumulate spikes per address
			sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) + delta_hi;
            sparse_sig = [sparse_sig; sig(1,TD.ts(spikeCount))];
        end
    else   % OFF spike
       plot([TD.ts(spikeCount)/1e6 TD.ts(spikeCount)/1e6], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),2)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),2)+0.4], 'k','Color','red');
        if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
            frame_data_off(TD.y(spikeCount)+1,TD.x(spikeCount)+1) = frame_data_off(TD.y(spikeCount)+1,TD.x(spikeCount)+1) + 1; % accumulate spikes per address          
			sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) - delta_lo;
            sparse_sig = [sparse_sig; sig(1,TD.ts(spikeCount))];
        end
    end 
end     
instrreset
%refreshdata
%ylim([0 191]);
%ylim([10 10]);
grid on
xlabel('Time (s)');
ylabel('Taxel Number');
%xytolin_customxy(TD.x,TD.y,12)

hold off;

subplot(3,1,2)
heatmap(frame_data_on);
colormap(jet(512))
title('ON spike count')

subplot(3,1,3)
heatmap(frame_data_off);
colormap(jet(512))   
title('OFF spike count')


% Set the width and height of the figure (in pixels)
width = 1200;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,['out\spike_frames_',name,'.fig'])
saveas(gcf,['out\spike_frames_',name,'.png'])
       
	   
%% Create a new figure for plotting sig and spline_sig together
figure;

% %for i = 1:mr_cnt
% for i=2:2
    i=1
    % Spline interpolation for each signal
    %spline_sig(i,:) = interp1(t, sig(i,:), t, 'linear');
    spline_sig(i,:) = interp1(TD.ts, sparse_sig, t, 'linear');
    for i = 1:size(spline_sig, 1)
    spline_sig(i, :) = fillmissing(spline_sig(i, :), 'linear', 'EndValues', 'nearest');
    end

    
%     % Create subplot for the signal and its spline interpolation
%     subplot(2, 2, i); % Adjust layout (2x2 grid for 4 subplots)
    
    % Plot the original signal
    plot(t, sig(i,:), 'LineWidth', 1.5); 
    hold on;
    
    % Plot the spline-interpolated signal
    plot(t, spline_sig(i,:), 'LineWidth', 1.5, 'LineStyle', '--'); 
    %plot(t, spline_sig_i, 'LineWidth', 1.5, 'LineStyle', '--'); 
    
    % Add labels and title for each subplot
    xlabel('Time (\mus)');
    ylabel('Amplitude');
    title(['Signal and Spline for Signal ', num2str(i)]);
    legend('Original Signal', 'Spline Interpolation');
    grid on;
    hold off;
% end

% Set the width and height of the figure (in pixels)
width = 1200;    % Specify your desired width here
height = 800;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

% Save the figure if needed
%saveas(gcf, ['out\sig_spline_', name, '.fig']);
saveas(gcf, ['out\sig_spline_', name, '.png']);



%% PLOT original reconstructed signal
figure;
% %for i=1:mr_cnt
% for i=2:2
    i=1
    %[psdEstimate_arr(i,:), freq_arr(i,:)] = pwelch(sig(i,:), windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided'); % psd not power matches better with sim		
    [psdEstimate_arr(i,:), freq_arr(i,:)] = pwelch(sig(i,:),  hann(windowSize), overlap, windowSize, samplingRate, 'psd', 'onesided'); % psd not power matches better with sim		
    
    %[psdEstimate_arr(i,:), freq_arr(i,:)] = pwelch(spline_sig(i,:), windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided'); % psd not power matches better with sim		

    %--added----
    % Frequency resolution (bin width)
    deltaF = samplingRate / windowSize;
    
    % Correct for windowing (Hann window coherent gain)
    windowGain = sum(hann(windowSize).^2) / windowSize; % Power correction factor

    psdEstimate = psdEstimate_arr(i,:) * deltaF / windowGain;

    % Compute Full Scale (FS) in dB
    FS_Vpp = 3; % Full scale in Vpp
    FS_RMS = FS_Vpp / (2 * sqrt(2)); % Convert Vpp to RMS
    FS_Power = FS_RMS^2; % Power corresponding to FS RMS
    FS_dB = 10 * log10(FS_Power); % Convert to dB

     psdEstimate = psdEstimate ./ FS_Power;

    %--added---

    %psdEstimate = psdEstimate_arr(i,:);
    freq = freq_arr(i,:);
    % Variables:
    % psdEstimate: Power Spectral Density (PSD) of the signal
    % freq: Frequency corresponding to the PSD

    % Find the index of the fundamental frequency component
    [~, fundIndex] = max(psdEstimate);
    fundIndex = find(freq == freqin);

    % Signal Power (at the fundamental frequency)
    signalPower = psdEstimate(fundIndex);

    % Define frequency bins for harmonics and noise
    harmonics = 10; % Number of harmonics to consider

    % Find harmonic frequencies and their indices
    harmonicIndices = zeros(harmonics,1);
    for h = 2:harmonics+1
        [~, harmonicIndices(h-1)] = min(abs(freq - h * freq(fundIndex)));
    end

    % Harmonic Power (sum of power at harmonic frequencies)
    harmonicPower = sum(psdEstimate(harmonicIndices));

    % Compute the noise power considering only frequencies greater than the fundamental frequency
    % Define the range for noise indices
    lowerFreqLimit = flo + 1/windowTime; %1 + 1/windowTime; %freqin * 1.2;
    upperFreqLimit = fbw; %lowerFreqLimit + fbw;
    
    % Find indices of frequencies between freqin*1.2 and freqin*1.2 + fbw
    noiseIndices = find(freq > lowerFreqLimit & freq <= upperFreqLimit);
    noiseIndices = setdiff(noiseIndices, harmonicIndices); % Exclude harmonic frequencies from noise
    noisePower = sum(psdEstimate(noiseIndices));

    % SNDR (Signal to Noise and Distortion Ratio)
    distortionPower = harmonicPower; % distortion from harmonics
    sndr = 10 * log10(signalPower / (noisePower + distortionPower));

    % SFDR (Spurious-Free Dynamic Range)
    spuriousPower = max(psdEstimate(noiseIndices)); % Find the highest spurious (non-signal) peak
    sfdr = 10 * log10(signalPower / spuriousPower);

    % SNR (Signal-to-Noise Ratio)
    snr = 10 * log10(signalPower / noisePower);

    enob = (sndr - 1.76)/6.02

    % Display the results
    fprintf('SNDR: %.2f dB\n', sndr);
    fprintf('SFDR: %.2f dB\n', sfdr);
    fprintf('SNR: %.2f dB\n', snr);
    fprintf('ENOB: %.2f bits\n', enob);

%     % Plot the Power Spectral Density (PSD)
%     subplot(2,2,i)
    plot(freq, 10*log10(psdEstimate), 'LineWidth', 1.5); % Convert PSD to dB scale for plotting
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (dBFS/Hz)');
    %title('Power Spectral Density (PSD)');
    set(gca,'XScale','log')
    grid on;
    xlim([0.5 1e3])

    % Annotate the plot with SNDR, SFDR, and SNR
    annotationX = 10; %min(freq) * 2; % x-position for text, adjust as needed
    annotationY = max(20*log10(psdEstimate)) - 10; % y-position for text, adjust as needed
    textStr = {
        sprintf('SNDR: %.2f dB', sndr), ...
        sprintf('SFDR: %.2f dB', sfdr), ...
        sprintf('SNR: %.2f dB', snr), ...
        sprintf('ENOB: %.2f bits', enob)
    };
    text(annotationX, annotationY, textStr, 'FontSize', 12, 'Color', 'red', 'VerticalAlignment', 'top');

    % Highlight the fundamental frequency and harmonic peaks
    hold on;
    plot(freq(fundIndex), 10*log10(psdEstimate(fundIndex)), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Highlight fundamental
    plot(freq(harmonicIndices), 10*log10(psdEstimate(harmonicIndices)), 'go', 'MarkerSize', 8, 'LineWidth', 2); % Highlight harmonics
    legend('PSD', 'Fundamental', 'Harmonics');
    hold off;
    
% end
   			   
micasplot

% Save the figure
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['for_paper\SPIKES\ADC\',name,'averaged_psd.fig']);
saveas(gcf, ['for_paper\SPIKES\ADC\',name,'averaged_psd.png']);

%% Function to fit polynomial and calculate error
% This function calculates the signal with given delta_hi and delta_lo, and returns the error
fit_error = @(deltas) calculate_fit_error(deltas, sig(2,:), ref_signal, t, TD, tduration, mr_cnt);
%fit_error = @(deltas) calculate_fit_error(deltas, spline_sig_i, ref_signal, t, TD, tduration, mr_cnt);


% Optimization of delta_hi and delta_lo with constraints
initial_guess = [delta_hi, delta_lo]; % Initial values of delta_hi and delta_lo
%==== CT_MODE = '1' (spikes reversed) =====%
lb = [-Inf, -Inf]; % Lower bounds for delta_hi and delta_lo (both must be >0)
%==== CT_MODE = '0'  =====%
%lb = [0, 0]; % Lower bounds for delta_hi and delta_lo (both must be >0)

%==== CT_MODE = '1' (spikes reversed) =====%
ub = [0, 0]; % Lower bounds for delta_hi and delta_lo (both must be >0)
%==== CT_MODE = '0'  =====%
%ub = [Inf, Inf]; % Upper bounds, can be set as needed

% Define optimization options
options = optimoptions('fmincon', 'Display', 'iter');

% Optimizing delta_hi and delta_lo with positivity constraint
opt_deltas = fmincon(@(deltas) calculate_fit_error(deltas, sig(2,:), ref_signal, t, TD, tduration, mr_cnt), ...
                     initial_guess, [], [], [], [], lb, ub, [], options);
%opt_deltas = fmincon(@(deltas) calculate_fit_error(deltas, spline_sig_i, ref_signal, t, TD, tduration, mr_cnt), ...
%                     initial_guess, [], [], [], [], lb, ub, [], options);

% Extract optimized values
delta_hi_opt = opt_deltas(1);
delta_lo_opt = opt_deltas(2);
%opt_deltas(1) = delta_hi; opt_deltas(2) = delta_lo;
% 
fprintf('Optimized delta_hi: %.5f V\n', delta_hi_opt);
fprintf('Optimized delta_lo: %.5f V\n', delta_lo_opt);


%opt_sig = spike_to_sig(opt_deltas,t,TD,tduration,mr_cnt);
opt_sig = spike_to_sig(opt_deltas,t,TD,tduration,mr_cnt); %to better fit p2p amplitude
%opt_sig = spike_to_sig_spline(opt_deltas,t,TD,tduration,mr_cnt); %to better fit p2p amplitude


%% PLOT optimized reconstructed signal
figure;
plot(t,opt_sig)

writematrix([t opt_sig],['out\',name,'-opt_sig.txt'])

figure;


[psdEstimate, freq] = pwelch(opt_sig, windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided'); % psd not power matches better with sim		

% Variables:
% psdEstimate: Power Spectral Density (PSD) of the signal
% freq: Frequency corresponding to the PSD

% Find the index of the fundamental frequency component
[~, fundIndex] = max(psdEstimate);
fundIndex = find(freq == freqin);

% Signal Power (at the fundamental frequency)
signalPower = psdEstimate(fundIndex);

% Define frequency bins for harmonics and noise
harmonics = 10; % Number of harmonics to consider

% Find harmonic frequencies and their indices
harmonicIndices = zeros(harmonics,1);
for h = 2:harmonics+1
    [~, harmonicIndices(h-1)] = min(abs(freq - h * freq(fundIndex)));
end

% Harmonic Power (sum of power at harmonic frequencies)
harmonicPower = sum(psdEstimate(harmonicIndices));

% Compute the noise power considering only frequencies greater than the fundamental frequency
% Define the range for noise indices
% lowerFreqLimit = freqin * 1.2;
% upperFreqLimit = lowerFreqLimit + fbw;
lowerFreqLimit = flo + 1/windowTime; %1 + 1/windowTime; %freqin * 1.2;
upperFreqLimit = fbw; %lowerFreqLimit + fbw;

% Find indices of frequencies between freqin*1.2 and freqin*1.2 + fbw
noiseIndices = find(freq > lowerFreqLimit & freq <= upperFreqLimit);
noiseIndices = setdiff(noiseIndices, harmonicIndices); % Exclude harmonic frequencies from noise
noisePower = sum(psdEstimate(noiseIndices));

% SNDR (Signal to Noise and Distortion Ratio)
distortionPower = harmonicPower; % distortion from harmonics
sndr = 10 * log10(signalPower / (noisePower + distortionPower));

% SFDR (Spurious-Free Dynamic Range)
spuriousPower = max(psdEstimate(noiseIndices)); % Find the highest spurious (non-signal) peak
sfdr = 10 * log10(signalPower / spuriousPower);

% SNR (Signal-to-Noise Ratio)
snr = 10 * log10(signalPower / noisePower);

enob = (sndr - 1.76)/6.02

% Display the results
fprintf('SNDR: %.2f dB\n', sndr);
fprintf('SFDR: %.2f dB\n', sfdr);
fprintf('SNR: %.2f dB\n', snr);
fprintf('ENOB: %.2f bits\n', enob);

% % Plot the Power Spectral Density (PSD)
% subplot(2,2,i)
plot(freq, 10*log10(psdEstimate), 'LineWidth', 1.5); % Convert PSD to dB scale for plotting
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density (PSD)');
set(gca,'XScale','log')
grid on;
xlim([0 10e3])

% Annotate the plot with SNDR, SFDR, and SNR
annotationX = 10; %min(freq) * 2; % x-position for text, adjust as needed
annotationY = max(20*log10(psdEstimate)) - 60; % y-position for text, adjust as needed
textStr = {
    sprintf('SNDR: %.2f dB', sndr), ...
    sprintf('SFDR: %.2f dB', sfdr), ...
    sprintf('SNR: %.2f dB', snr), ...
    sprintf('ENOB: %.2f bits', enob)
};
text(annotationX, annotationY, textStr, 'FontSize', 12, 'Color', 'red', 'VerticalAlignment', 'top');

% Highlight the fundamental frequency and harmonic peaks
hold on;
plot(freq(fundIndex), 10*log10(psdEstimate(fundIndex)), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Highlight fundamental
plot(freq(harmonicIndices), 10*log10(psdEstimate(harmonicIndices)), 'go', 'MarkerSize', 8, 'LineWidth', 2); % Highlight harmonics
legend('PSD', 'Fundamental', 'Harmonics');
hold off;


micasplot

% Save the figure
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['for_paper\SPIKES\ADC\',name,'averaged_psd_opt.fig']);
saveas(gcf, ['for_paper\SPIKES\ADC\',name,'averaged_psd_opt.png']);


function sig = spike_to_sig(deltas,t,TD,tduration,mr_cnt)

    delta_hi = deltas(1);
    delta_lo = deltas(2);

    sig = zeros(1,length(t)); 
    idx_within_time = find(TD.ts<tduration*1e6);
     for spikeCount = 1:max(idx_within_time)
        if TD.p(spikeCount)==1  % ON spike
           if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
                sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) + delta_hi;
           end
        else   % OFF spike
            if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
                sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) - delta_lo;
            end
        end 
     end  


end

function spline_sig_i = spike_to_sig_spline(deltas,t,TD,tduration,mr_cnt)

    delta_hi = deltas(1);
    delta_lo = deltas(2);

    sig = zeros(1,length(t)); 
    sparse_sig = [0];
    idx_within_time = find(TD.ts<tduration*1e6);
     for spikeCount = 1:max(idx_within_time)
        if TD.p(spikeCount)==1  % ON spike
           if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
                sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) + delta_hi;
                sparse_sig = [sparse_sig; sig(1,TD.ts(spikeCount))];
           end
        else   % OFF spike
            if TD.ts(spikeCount) > 0 %1000000 % 1 second; ignore initial spikes
                sig(1,TD.ts(spikeCount):length(t)) = sig(1,TD.ts(spikeCount):length(t)) - delta_lo;
                sparse_sig = [sparse_sig; sig(1,TD.ts(spikeCount))];
            end
        end 
    end  
     spline_sig_i = interp1(TD.ts, sparse_sig, t, 'spline');

end

function error = calculate_fit_error_tnorm(deltas, sig, ref_signal, t, TD, tduration, mr_cnt)
    delta_hi = deltas(1);
    delta_lo = deltas(2);
    
    % Example of updating the signal with delta_hi and delta_lo
    % Adjust this logic based on how delta_hi and delta_lo are used in your signal
    % Here, we assume that delta_hi and delta_lo modulate the signal's spikes
    %modulated_sig = sig; % Start with the original signal and modify
    modulated_sig = spike_to_sig(deltas,t,TD, tduration, mr_cnt);
    % Polynomial fitting (replace with your actual fitting method)
    
    % Centering and scaling the x-data
    t_mean = mean(t);
    t_std = std(t);
    t_scaled = (t - t_mean) / t_std; % Center and scale time data
    
    % Fit a lower-order polynomial, e.g., 2nd or 3rd degree
    p = polyfit(t_scaled, modulated_sig, 5); % Fit a 2nd-order polynomial
    
    % Generate the fitted polynomial signal
    fitted_signal = polyval(p, t_scaled);
    
    % Calculate the error between the fitted signal and reference sinusoidal signal
    error = sum((fitted_signal - ref_signal).^2) % Sum of squared differences
end

function error = calculate_fit_error(deltas, sig, ref_signal, t, TD, tduration, mr_cnt)
    delta_hi = deltas(1);
    delta_lo = deltas(2);
    
    % Example of updating the signal with delta_hi and delta_lo
    % Adjust this logic based on how delta_hi and delta_lo are used in your signal
    % Here, we assume that delta_hi and delta_lo modulate the signal's spikes
    %modulated_sig = sig; % Start with the original signal and modify
    modulated_sig = spike_to_sig(deltas,t,TD, tduration, mr_cnt);
    % Polynomial fitting (replace with your actual fitting method)
    
%     % Centering and scaling the x-data
%     t_mean = mean(t);
%     t_std = std(t);
%     t_scaled = (t - t_mean) / t_std; % Center and scale time data
    
    % Fit a lower-order polynomial, e.g., 2nd or 3rd degree
    p = polyfit(t, modulated_sig, 5); % Fit a 2nd-order polynomial
    
    % Generate the fitted polynomial signal
    fitted_signal = polyval(p, t);
    
    % Calculate the error between the fitted signal and reference sinusoidal signal
    error = sum((fitted_signal - ref_signal).^2) % Sum of squared differences
end

