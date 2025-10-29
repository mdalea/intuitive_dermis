close all

% Define input values
no_of_channels = [16, 70, 70, 19, 10, 4, 4];
adc_bits = [24, 12, 12, 12, 14, 6, 6];
sampling_rate = [300, 30, 30, 4.4e3, 390, 2*50, 2*50];
%power = 1e6 * [42e-3, 5*10e-3, 5*10e-3, 5*5e-3, 5*258e-6, 7.6e-6*4+0.6e-3*3+2e-3*6, 7.6e-6*4]; % Power in µW
%power = 1e6 * [42e-3, 5*10e-3, 5*10e-3, 5*5e-3, 5*258e-6+3.3*100e-3, 7.6e-6*4+0.6e-3*3+2e-3*6, 7.6e-6*4]; % Power in µW
power = 1e6 * [42e-3, 5*10e-3, 5*10e-3, 5*5e-3, 5*258e-6+3.3*100e-3, 10.6e-6*4+0.6e-3*3+2e-3*6, 10.6e-6*4]; % Power in µW
%labels = {sprintf('[2]''20\nADS1258'), sprintf('[14]''14\nDSA9210'),  sprintf('[15]''15\nDSA9210'), ...
%          sprintf('[16]''15\nPIC\muC'), sprintf('[17]''23\nADPD1080'), ...
%          sprintf('This Work\na-IGZO TFT'), sprintf('This Work\na-IGZO TFT\n(w/o AER Power Consumption)')};
labels = {sprintf('[3]\nADS1258'), sprintf('[4]\nDSA9210'),  sprintf('[5]\nDSA9210'), ...
          sprintf('[16]\nPIC\\muC'), sprintf('[17]\nADPD1080\n+ Teensy 4.1'), ...
          sprintf('This Work\na-IGZO TFT'), sprintf('This Work\na-IGZO TFT\n(w/o AER Power\nConsumption)')};
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
    if i == 1
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');    
    elseif i == 2
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    elseif i == 3
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');        
    elseif i == 5    
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');        
    elseif i == 6    
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');
    elseif i == 7    
        text(log10(raw_data_rate(i)), log10(power(i)), labels{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'r');    
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
ylabel('Readout Power Consumption (\muW)');

% Fix x-axis ticks manually for log scale
xtick_values = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]; % Define tick values for x-axis
xticks(log10(xtick_values)); % Set logarithmic tick positions for x-axis
xticklabels(compose('%.0e', xtick_values)); % Format labels for x-axis in scientific notation

% Create the legend and increase the icon size by 3x
lgd = legend('Envelope Detection', 'Transforms', 'Cross-Correlation', 'Classifier', 'Friction Coefficient', ...
             sprintf('Friction Coefficient\n(Grasp State Aware)'), sprintf('Friction Coefficient\n(Grasp State Aware)'), 'Location', 'southeast');
lgd.FontSize = 12;  % Adjust the font size if needed
legendItems = lgd.EntryContainer.Children;  % Get the legend items
for i = 1:length(legendItems)
    legendItems(i).MarkerSize = 18;  % Increase the marker size (3x the default size)
end

% Set axis limits to ensure text fits
micasplot
% Set y-axis limits to 100*max(power)
ylim([0, log10(100 * max(power))]);  % Set the y-axis limit to 100 times the maximum power value

% Save the ZOH plot
width = 800;
height = 500;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['./power_vs_datarate.png']);
