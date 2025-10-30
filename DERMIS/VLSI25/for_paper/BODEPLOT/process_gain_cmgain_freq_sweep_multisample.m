% READ output of keysight_rigol_gain_freqsweep.py over N iterations for averaging
close all;
clear all;

% Define run IDs and corresponding run_iter for each sample
samples_gain = {
    'sample_3_vtune_0p2_vcmin_1p513_vcm_1p5_vddafe_3p0_ibias_20n_fixed', 
    'sample_4_ibias_50n'
};
run_iters_gain = [1, 1]; % Number of iterations for each sample

% Define colors for the samples
colors = {'r', 'b'}; % Red for Sample 1, Blue for Sample 2

figure;

% Array to store handles for legend
gain_handles = [];

% Loop through each sample for Gain data
for sample_idx = 1:length(samples_gain)
    run_id = samples_gain{sample_idx};
    baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-'];
    run_iter = run_iters_gain(sample_idx); % Get the corresponding run_iter

    freq = []; gain = [];
    
    % Read and process data for the current sample
    for i = 1:run_iter
        csvFileName = [baseName, num2str(i-1), '.csv'];
        csv_val = csvread(csvFileName, 1, 0); % Skip the first line of column names
        
        freq(i, :) = csv_val(:, 1);
        gain(i, :) = csv_val(:, 4);
    end

    freq = freq';
    gain = gain';

    % Compute the mean across iterations
    mean_freq = mean(freq, 2);
    mean_gain = mean(gain, 2);

    % Plot the gain for the current sample
    subplot(2, 1, 1);
    hold on;
    % Add line plot with color and save handle
    h = plot(mean_freq, 20 * log10(mean_gain), '-', 'Color', colors{sample_idx}, 'LineWidth', 1.5);
    gain_handles = [gain_handles, h]; % Store handle for the legend
    % Add scatter plot with squares and matching color
    scatter(mean_freq, 20 * log10(mean_gain), 100, 's', 'filled', 'MarkerFaceColor', colors{sample_idx});
    
    % Set axis properties
    set(gca, 'XScale', 'log');
    grid on;
    ylim([-30 40]);
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
end

% Finalize the first subplot with legend
legend(gain_handles, {'Sample 1', 'Sample 2'}, 'location', 'best');
title('Gain vs Frequency');
hold off;

%% Repeat for CM Gain data
samples_cmgain = {
    'sample_3_vtune_0p2_vcmin_1p513_vcm_1p5_vddafe_3p0_ibias_40n_fixed_2', 
    'sample_4_ibias_50n'
};
run_iters_cmgain = [3, 3]; % Number of iterations for each sample in CM gain

% Array to store handles for CM Gain legend
cmgain_handles = [];

for sample_idx = 1:length(samples_cmgain)
    run_id = samples_cmgain{sample_idx};
    baseName = ['output\ID_', run_id, '_cmgain__gain_vs_fin_sweep-'];
    run_iter = run_iters_cmgain(sample_idx); % Get the corresponding run_iter

    freq_cm = []; gain_cm = [];
    
    % Read and process CM gain data for the current sample
    for i = 1:run_iter
        csvFileName = [baseName, num2str(i-1), '.csv'];
        csv_val = csvread(csvFileName, 1, 0); % Skip the first line of column names
        
        freq_cm(i, :) = csv_val(:, 1);
        gain_cm(i, :) = csv_val(:, 4);
    end

    freq_cm = freq_cm';
    gain_cm = gain_cm';

    % Compute the mean across iterations
    mean_freq_cm = mean(freq_cm, 2);
    mean_gain_cm = mean(gain_cm, 2);

    % Plot the CM gain for the current sample
    subplot(2, 1, 2);
    hold on;
    % Add line plot with color and save handle
    h = plot(mean_freq_cm, 20 * log10(mean_gain_cm), '-', 'Color', colors{sample_idx}, 'LineWidth', 1.5);
    cmgain_handles = [cmgain_handles, h]; % Store handle for the legend
    % Add scatter plot with squares and matching color
    scatter(mean_freq_cm, 20 * log10(mean_gain_cm), 100, 's', 'filled', 'MarkerFaceColor', colors{sample_idx});

    % Set axis properties
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('CM Gain (dB)');
end

% Finalize the second subplot with legend
legend(cmgain_handles, {'Sample 1', 'Sample 2'}, 'location', 'best');
title('CM Gain vs Frequency');
hold off;

micasplot
% Save the figure
% Set the width and height of the figure (in pixels)
width = 600; %800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, 'output\superimposed_gain_cmgain.fig');
saveas(gcf, 'output\superimposed_gain_cmgain.png');
