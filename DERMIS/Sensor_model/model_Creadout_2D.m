clear all
close all
clc

%% ==== SENSOR/READOUT PARAMS =====
N   = 5;
Cf  = 10e-15;
VDD = 3.3;
Vocm = VDD/2; % source of nonlinearity for VCM measurements

%% === NO SHEAR FORCE; NORMAL FORCE ONLY
x0 = 50e-6;   % overlap at shear force = 0  (m)
z0 = 125e-6;  % thickness at normal force = 0 (m)
y0 = 125e-6;  % constant height overlap (m)

deltaZ_max = 100e-6; 
deltaZ_res = deltaZ_max./2^N./2; 
deltaZ      = 0:deltaZ_res:deltaZ_max; % full sweep of deltaZ

deltaX_max = 50e-6;
deltaX_res = deltaX_max./2^N;
deltaX     = -deltaX_max:deltaX_res:deltaX_max;

% Create grid
[DeltaX, DeltaZ] = meshgrid(deltaX, deltaZ);

% Dielectric constants
eps0      = 8.85e-12;     % F/m
eps_pdms  = 2.77;         
eps_eff   = eps_pdms * eps0; % effective permittivity (F/m)

%% ==== Capacitances ====
% Csp (formerly Cs1)
Cspa = (eps_eff .* (x0 - DeltaX) .* y0) ./ (z0 - DeltaZ);
Cspb = (eps_eff .* (x0 - DeltaX) .* y0) ./ (z0 - DeltaZ);
Csp  = (Cspa .* Cspb) ./ (Cspa + Cspb);
Csp(isnan(Csp)) = 0;  % remove 0/0 singularities

% Csn (formerly Cs2)
Csna = (eps_eff .* (x0 + DeltaX) .* y0) ./ (z0 - DeltaZ);
Csnb = (eps_eff .* (x0 + DeltaX) .* y0) ./ (z0 - DeltaZ);
Csn  = (Csna .* Csnb) ./ (Csna + Csnb);
Csn(isnan(Csn)) = 0;

%% ==== 2D PLOTS: Csp and Csn vs Delta-X for selected Delta-Z ==== 
figure('Color','w');
hold on;

% Pick 5 evenly spaced indices from deltaZ range
num_steps = 5;
z_indices = round(linspace(1, length(deltaZ), num_steps));

colors = lines(num_steps); % distinct colors for each ΔZ slice

% Plot Csp first (solid lines)
hCsp = gobjects(num_steps,1);
for i = 1:num_steps
    idx = z_indices(i);
    hCsp(i) = plot(deltaX*1e6, Csp(idx,:)*1e15, ...
        'Color', colors(i,:), 'LineStyle','-', 'LineWidth',1.6); % fF and µm
end

% Then plot Csn (dashed lines)
hCsn = gobjects(num_steps,1);
for i = 1:num_steps
    idx = z_indices(i);
    hCsn(i) = plot(deltaX*1e6, Csn(idx,:)*1e15, ...
        'Color', colors(i,:), 'LineStyle','--', 'LineWidth',1.6); % fF and µm
end

grid on
xlabel('{\Delta}X (µm)')
ylabel('Capacitance (fF)')
title('C_{sp} (solid) and C_{sn} (dashed) vs {\Delta}X at selected {\Delta}Z')

% Build legend: all Csp entries first, then Csn
legCsp = arrayfun(@(z) sprintf('C_{sp}, \\DeltaZ = %.1f µm', deltaZ(z)*1e6), z_indices, 'UniformOutput', false);
legCsn = arrayfun(@(z) sprintf('C_{sn}, \\DeltaZ = %.1f µm', deltaZ(z)*1e6), z_indices, 'UniformOutput', false);
hLeg = legend([hCsp; hCsn], [legCsp, legCsn], 'Location','bestoutside');

% --- Axis arrow annotations (red) ---
% X-axis arrow for Shear Force (↔)
annotation('textarrow',[0.25 0.75],[0.12 0.12], ...
    'String','Shear Force','HeadStyle','plain','LineWidth',1.2, 'Color','r', 'TextColor','r');
annotation('textarrow',[0.75 0.25],[0.12 0.12], ...
    'String','','HeadStyle','plain','LineWidth',1.2, 'Color','r');

% Y-axis arrow for Normal Force (↕)
annotation('textarrow',[0.08 0.08],[0.25 0.75], ...
    'String','Normal Force','HeadStyle','plain','LineWidth',1.2, 'Color','r', 'TextColor','r');
annotation('textarrow',[0.08 0.08],[0.75 0.25], ...
    'String','','HeadStyle','plain','LineWidth',1.2, 'Color','r');

% --- Parameter text below legend ---
param_str = sprintf(['x_0 = %.1f µm\n' ...
                     'y_0 = %.1f µm\n' ...
                     'z_0 = %.1f µm\n' ...
                     '\\epsilon_{PDMS}{\\cdot}\\epsilon_0 = %.3g F/m'], ...
                     x0*1e6, y0*1e6, z0*1e6, eps_eff);

% Place annotation box just below the legend position
legPos = get(hLeg,'Position'); % [x y w h]
annotation('textbox',[legPos(1) legPos(2)-0.15 legPos(3) 0.12], ...
    'String', param_str, 'FitBoxToText','on', ...
    'BackgroundColor',[1 1 1], 'EdgeColor',[0.7 0.7 0.7], ...
    'FontSize',9, 'LineWidth',1);


micasplot

% Save the  plot
width = 800;
height = 400;
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf, ['./Csensors_DR.png']);