close all

% Define input values
latency = 1e3 * [40e-3, 64e-3, 1/20, 10e-3, 1/40, 10e-3]; % in ms
spatialres = [2.36, 3.4, 3.4, 2, 7.5, 0.6];

labels = {sprintf('[3]\n(Time Window)'), sprintf('[4]'),  sprintf('[5]'), ...
          sprintf('[16]\n(Time Window)'), sprintf('[17]'), ...
          sprintf('This Work\n(Spike Counting Period)')};
markers = {'p', 's', 'o', 'd', '^', 'v'}; % Marker types

% Plot the results as a scatter plot
figure;
hold on;
for i = 1:length(latency)
    marker_size = 500;
    scatter(spatialres(i), latency(i), marker_size, 'filled', markers{i});

    if i == 4
        text(spatialres(i), latency(i), labels{i}, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'FontSize', 10);
    elseif i==6
        text(spatialres(i), latency(i), labels{i}, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
            'FontSize', 10, 'Color', 'r');
      
    else
        text(spatialres(i), latency(i), labels{i}, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
            'FontSize', 10);
    end
end

% Add dashed horizontal line at 60 ms with label
yline(60, '--', sprintf('Total Latency\n(Humans)'), ...
      'FontSize', 15, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');

hold off;

% Set axis limits to ensure text fits
xlim([min(spatialres) - 0.5, max(spatialres) + 5]);
ylim([-10, max(latency) + 20]);

grid on;
xlabel('Spatial Resolution (mm)');
ylabel('Slip Detection Latency (ms)');

% Create the legend and increase the icon size by 3x
lgd = legend('Envelope Detection', 'Transforms', 'Cross-Correlation', 'Classifier', 'Friction Coefficient', ...
             sprintf('Friction Coefficient\n(Grasp State Aware)'), 'Location', 'northeast');
lgd.FontSize = 12;  % Adjust the font size if needed
legendItems = lgd.EntryContainer.Children;  % Get the legend items
for i = 1:length(legendItems)
    legendItems(i).MarkerSize = 18;  % Increase the marker size (3x the default size)
end

micasplot

% Save the ZOH plot
width = 800;
height = 400;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['./latency_vs_spatialres.png']);


