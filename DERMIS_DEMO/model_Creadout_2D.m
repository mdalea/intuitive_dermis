clear all
close all
clc

%% ==== SENSOR/READOUT PARAMS =====
N   = 5;
Cf  = 10e-15;
VDD = 3.3;
Vocm = VDD/2;

%% === NO SHEAR FORCE; NORMAL FORCE ONLY
x0 = 50e-6;
z0 = 125e-6;
y0 = 125e-6;

deltaZ_max = 100e-6; 
deltaZ_res = deltaZ_max./2^N./2; 
deltaZ      = 0:deltaZ_res:deltaZ_max;

deltaX_max = 50e-6;
deltaX_res = deltaX_max./2^N;
deltaX     = -deltaX_max:deltaX_res:deltaX_max;

[DeltaX, DeltaZ] = meshgrid(deltaX, deltaZ);

eps0      = 8.85e-12;     % F/m
eps_pdms  = 2.77;         
eps_eff   = eps_pdms * eps0;

%% ==== Capacitances ====
Cspa = (eps_eff .* (x0 - DeltaX) .* y0) ./ (z0 - DeltaZ);
Cspb = (eps_eff .* (x0 - DeltaX) .* y0) ./ (z0 - DeltaZ);
Csp  = (Cspa .* Cspb) ./ (Cspa + Cspb);
Csp(isnan(Csp)) = 0;

Csna = (eps_eff .* (x0 + DeltaX) .* y0) ./ (z0 - DeltaZ);
Csnb = (eps_eff .* (x0 + DeltaX) .* y0) ./ (z0 - DeltaZ);
Csn  = (Csna .* Csnb) ./ (Csna + Csnb);
Csn(isnan(Csn)) = 0;

%% ==== Plot 1: Csp and Csn vs Delta-X for selected Delta-Z ==== 
figure('Color','w');
hold on;

num_steps = 5;
z_indices = round(linspace(1, length(deltaZ), num_steps));
colors = lines(num_steps);

hPlots = [];
legText = {};

for i = 1:num_steps
    idx = z_indices(i);
    h1 = plot(deltaX*1e6, Csp(idx,:)*1e15, ...
        'Color', colors(i,:), 'LineStyle','-', 'LineWidth',1.6);
    h2 = plot(deltaX*1e6, Csn(idx,:)*1e15, ...
        'Color', colors(i,:), 'LineStyle','--', 'LineWidth',1.6);
    hPlots = [hPlots, h1, h2];
    legText = [legText, {sprintf('ΔZ = %.1f µm, Csp', deltaZ(idx)*1e6)}, ...
                         {sprintf('ΔZ = %.1f µm, Csn', deltaZ(idx)*1e6)}];
end

grid on
xlabel('{\Delta}X (µm)')
ylabel('Capacitance (fF)')
legend(hPlots, legText, 'Location','bestoutside');

micasplot
set(gcf, 'Position', [100, 100, 800, 400]);
saveas(gcf, './Csensors_DR.png');

%% ==== Plot 2: (Csp + Csn) vs Delta-X for the same Delta-Z slices ====
figure('Color','w');
hold on;

hSum = gobjects(num_steps,1);
sumLeg = cell(1,num_steps);

for i = 1:num_steps
    idx = z_indices(i);
    Csum = (Csp(idx,:) + Csn(idx,:)) * 1e15; % fF
    hSum(i) = plot(deltaX*1e6, Csum, ...
        'Color', colors(i,:), 'LineStyle','-', 'LineWidth',1.8);
    sumLeg{i} = sprintf('ΔZ = %.1f µm, Csp + Csn', deltaZ(idx)*1e6);
end

grid on
xlabel('{\Delta}X (µm)')
ylabel('Capacitance (fF)')
legend(hSum, sumLeg, 'Location','bestoutside');

micasplot
set(gcf, 'Position', [120, 520, 800, 400]);
saveas(gcf, './Csensors_sum_DR.png');

