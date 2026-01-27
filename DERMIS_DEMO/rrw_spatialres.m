close all

% Input values (latency removed, only spatial resolution kept)
spatialres = [2.36, 3.4, 3.4, 2, 7.5, 0.6];

labels = {sprintf('[4]'), sprintf('[5]'),  sprintf('[6]'), ...
          sprintf('[3]'), sprintf('[18]'), ...
          sprintf('This\nWork\n(DERMIS)')};

% Positions along x-axis for each square
N = length(spatialres);
xpos = 1:N;
ypos = zeros(1, N);  % all on the same horizontal line

figure;
hold on;

% Scale factor for marker size (tune if needed)
base_size = 4000;

for i = 1:N
    % Marker size proportional to spatial resolution (area ∝ spatialres)
    marker_size = base_size * (spatialres(i) / max(spatialres))^2

    % Plot filled square
    scatter(xpos(i), ypos(i), marker_size, 's', 'filled');

    % Put label on/near the square
    text(xpos(i), ypos(i), labels{i}, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
        'FontSize', 10);
end

hold off;

% Make plot look nicer
xlim([0.5, N + 0.5]);
ylim([-1, 2]);                % Just to give some room for the labels
set(gca, 'YTick', []);        % No meaningful y-axis
set(gca, 'XTick', xpos);      % Optional: ticks under each square
set(gca, 'XTickLabel', []);   % Hide x tick labels since we use text()

grid on;
xlabel('');
ylabel('');  % No latency axis anymore
title('Spatial Resolution Comparison (Square Area ∝ Spatial Resolution)');

micasplot

% Save the plot
width = 800;
height = 400;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, './spatialres_squares.png');

