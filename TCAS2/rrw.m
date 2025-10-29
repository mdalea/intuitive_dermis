close all

% Define input values
no_of_channels = [16, 70, 70, 19, 10, 4, 4];
adc_bits = [24, 12, 12, 12, 14, 6, 6];
sampling_rate = [300, 30, 30, 4.4e3, 390, 2*50, 2*50];
power = 1e6 * [42e-3, 5*10e-3, 5*10e-3, 5*5e-3, 5*258e-6, 7.6e-6*4+0.6e-3*3+2e-3*6, 7.6e-6*4]; % Power in µW
labels = {sprintf('[2]''20\nADS1258'), sprintf('[14]''14\nDSA9210'),  sprintf('[15]''15\nDSA9210'), ...
          sprintf('[16]''15\nPIC\muC'), sprintf('[17]''23\nADPD1080'), ...
          sprintf('This Work\na-IGZO TFT'), sprintf('This Work\na-IGZO TFT\n(w/o AER Power Consumption)')};
markers = {'p', 's', 'o', 'd', '^', 'v', 'v'}; % Marker types

% Compute raw data rate
raw_data_rate = no_of_channels .* adc_bits .* sampling_rate;

% Plot the results as a scatter plot
figure;
hold on;
for i = 1:length(raw_data_rate)
    marker_size = 500;
    if i == 3
        marker_size = 100; % Smaller size for 3rd item
    end
    scatter(log10(raw_data_rate(i)), log10(power(i)), marker_size, 'filled', markers{i});
    
    if i == 2
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    else
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end
hold off;

% Set y-axis to display actual power values in scientific notation
ytick_values = [1 10 100 1e3 1e4 1e5 1e6]; % Define tick values
yticks(log10(ytick_values)); % Set logarithmic tick positions
yticklabels(compose('%.0e', ytick_values)); % Format labels in scientific notation

grid on;
xlabel('Raw Data Rate (bps)');
ylabel('Power Consumption (\muW)');
title('System Power Consumption vs. Raw Data Rate');
legend('Envelope Detection', 'Transforms', 'Cross-Correlation', 'Classifier', 'Friction Coefficient', ...
       'Friction Coefficient (on-chip)', 'Friction Coefficient (on-chip)', 'Location', 'southeast');


micasplot