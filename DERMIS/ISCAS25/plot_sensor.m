% Clear workspace and command window
clear;
clc;

% Specify the filename format
% You can change the filename here to the one you want to analyze
filename = 'arduino_yy_YshearForce1_F1.3_TR167.csv';
%filename = 'arduino_xx_XshearForce1_F1.3_TR171.csv';

% Read the CSV file (assuming no header in the file)
data = readmatrix(filename);

% Extract capacitance values (ignore time)
Csp = data(:, 2);     % Second column is Csp (capacitance)

% Separate Csp into Csp1 and Csp2 (alternating entries)
Csp1 = Csp(1:2:end);  % Odd indexed values
Csp2 = Csp(2:2:end);  % Even indexed values

% Filter: Retain only values where 0 < Csp1 < 1 and 0 < Csp2 < 1
Csp1 = Csp1(Csp1 > 0 & Csp1 < 1);
Csp2 = Csp2(Csp2 > 0 & Csp2 < 1);

% Create time axis for Csp1 and Csp2 (in seconds, each index = 11ms)
time_step = 0.011;  % Time step in seconds (11 ms)
time_Csp1 = (0:length(Csp1)-1) * time_step;  % Time for Csp1
time_Csp2 = (0:length(Csp2)-1) * time_step;  % Time for Csp2

% ---- Apply Low Pass Filter (LPF) to Csp1 and Csp2 ----

% Sampling rate (assuming the original sampling frequency is 10 kHz)
fs = 10000;  % Sampling frequency in Hz
fc = 1000;   % Cutoff frequency (1 kHz)

% Design a 2nd-order Butterworth low-pass filter with 1 kHz cutoff
[b, a] = butter(2, fc/(fs/2));  % Normalized cutoff frequency

% Apply the filter to Csp1 and Csp2
Csp1_filt = filtfilt(b, a, Csp1);  % Zero-phase filtering to avoid phase distortion
Csp2_filt = filtfilt(b, a, Csp2);

% ---- Plot the filtered data and calculations ----

% Create new figure for filtered data plots
figure;

% Subplot 1: Csp1_filt and Csp2_filt on the same plot (Time axis in seconds)
subplot(4, 1, 1);
plot(time_Csp1, Csp1_filt, 'r-', 'DisplayName', 'Csp'); hold on;
plot(time_Csp2, Csp2_filt, 'b--', 'DisplayName', 'Csn');
%xlabel('Time (s)');
%ylabel('(pF)');
legend('show');
title('Raw Csp and Csn (pF)');
xlim([10 40])
grid on

% Subplot 2: Csp1_filt - Csp2_filt
subplot(4, 1, 2);
minLength = min(length(Csp1_filt), length(Csp2_filt)); % To handle cases where lengths differ
plot(time_Csp1(1:minLength), Csp1_filt(1:minLength) - Csp2_filt(1:minLength), 'g-');
%%xlabel('Time (s)');
%ylabel('Csp - Csn');
title('Raw differential capacitance: Csp - Csn (pF)');
xlim([10 40])
grid on

% Subplot 3: Csp1_filt + Csp2_filt
subplot(4, 1, 3);
plot(time_Csp1(1:minLength), (Csp1_filt(1:minLength) + Csp2_filt(1:minLength))/0.3, 'm-');
%xlabel('Time (s)');
ylabel(['Normal\newlineForce\newline[N]']);
title('Csp + Csn - Sensitivity: 0.3pF/N');
xlim([10 40])
grid on

% Subplot 4: (Csp1_filt - Csp2_filt) / (Csp1_filt + Csp2_filt)
subplot(4, 1, 4);
ratio_filt = (Csp1_filt(1:minLength) - Csp2_filt(1:minLength)) ./ (Csp1_filt(1:minLength) + Csp2_filt(1:minLength));
plot(time_Csp1(1:minLength), ratio_filt, 'k-');
xlabel('Time (s)');
ylabel('Shear\newlineForce\newline[pF/pF]');
title('Shear Force: (Csp - Csn) / (Csp + Csn)');


% ---- Add New Subplot with 1 Hz LPF ----

% Design a 2nd-order Butterworth low-pass filter with 1 Hz cutoff
fc_1Hz = 100;  % Cutoff frequency (1 Hz)
[b_1Hz, a_1Hz] = butter(2, fc_1Hz/(fs/2));  % Normalized cutoff frequency

% Apply the 1 Hz filter to Csp1_filt and Csp2_filt
Csp1_filt_1Hz = filtfilt(b_1Hz, a_1Hz, Csp1_filt);  % Zero-phase filtering to avoid phase distortion
Csp2_filt_1Hz = filtfilt(b_1Hz, a_1Hz, Csp2_filt);

% % Subplot 5: Ratio with 1 Hz cutoff filter applied
% subplot(5, 1, 5);
hold on;
ratio_filt_1Hz = (Csp1_filt_1Hz(1:minLength) - Csp2_filt_1Hz(1:minLength)) ./ (Csp1_filt_1Hz(1:minLength) + Csp2_filt_1Hz(1:minLength));
plot(time_Csp1(1:minLength), ratio_filt_1Hz, 'b-');

xlim([10 40])
grid on

micasplot


% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,['Sensor prototype\',filename,'-capvals','.fig'])
saveas(gcf,['Sensor prototype\',filename,'-capvals','.png'])


% xlabel('Time (s)');
% ylabel('Ratio (1 Hz LPF)');
% title('Ratio with 1 Hz LPF');
% 
% % Adjust the figure layout for clarity
% sgtitle(['Filtered Capacitance Analysis with <1kHz and <1Hz LPF for ', filename]);
% 
% % Optionally, add a text annotation at the bottom of the figure
% dim = [0.13, 0.01, 0.8, 0.1];  % [x, y, width, height] in normalized figure units
% annotation('textbox', dim, 'String', ['Analyzed CSV File: ', filename], 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'center');
