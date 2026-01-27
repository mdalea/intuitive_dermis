clear all
close all
clc

%% ==== SENSOR/READOUT PARAMS =====
N   = 5;
Cf  = 10e-15;     % F
VDD = 3.3;        % V
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

%% ==== Vout(diff) ====
% Vout_diff = (VDD/Cf) * (Csp - Csn)
Vout_diff = (VDD./Cf) .* (Csp - Csn);   % V

%% ==== Print equivalent (Csp - Csn) for Vout(LSB/2) = 45 mV ====
Vout_LSB2 = 45e-3;                      % V
dC_eq     = (Vout_LSB2 .* Cf) ./ VDD;   % F

fprintf('\n==== Equivalent capacitance difference for Vout(LSB/2) = 45 mV ====\n');
fprintf('Cf = %.3f fF, VDD = %.3f V\n', Cf*1e15, VDD);
fprintf('Csp - Csn = %.6f fF  (= %.3f aF)\n\n', dC_eq*1e15, dC_eq*1e18);

%% ==== Extra dynamic range measurement: 20log10(max/min) ====
Vpp_max = 6.0;        % Vpp differential
Vpp_min = 45e-3;      % Vpp differential

DR_dB = 20*log10(Vpp_max./Vpp_min);

fprintf('==== Dynamic Range (given) ====\n');
fprintf('max = %.3f Vpp diff, min = %.3f Vpp diff\n', Vpp_max, Vpp_min);
fprintf('DR = 20log10(max/min) = %.3f dB\n\n', DR_dB);

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

%% ==== Plot 3: Vout(diff) vs Delta-X for the same Delta-Z slices (same colors) ====
figure('Color','w');
hold on;

hVx = gobjects(num_steps,1);
legVx = cell(1,num_steps);

for i = 1:num_steps
    idx = z_indices(i);
    hVx(i) = plot(deltaX*1e6, Vout_diff(idx,:), ...
        'Color', colors(i,:), 'LineStyle','-', 'LineWidth',1.8);
    legVx{i} = sprintf('ΔZ = %.1f µm', deltaZ(idx)*1e6);
end

grid on
xlabel('{\Delta}X (µm)')
ylabel('V_{out,diff} (V)')
legend(hVx, legVx, 'Location','bestoutside');

micasplot
set(gcf, 'Position', [140, 940, 800, 400]);
saveas(gcf, './Voutdiff_vs_DeltaX.png');

%% ==== Plot 4: Vout(diff) vs Delta-Z for selected Delta-X slices (same colors) ====
figure('Color','w');
hold on;

x_indices = round(linspace(1, length(deltaX), num_steps));

hVz = gobjects(num_steps,1);
legVz = cell(1,num_steps);

for i = 1:num_steps
    idx = x_indices(i);
    hVz(i) = plot(deltaZ*1e6, Vout_diff(:,idx), ...
        'Color', colors(i,:), 'LineStyle','-', 'LineWidth',1.8);
    legVz{i} = sprintf('ΔX = %.1f µm', deltaX(idx)*1e6);
end

grid on
xlabel('{\Delta}Z (µm)')
ylabel('V_{out,diff} (V)')
legend(hVz, legVz, 'Location','bestoutside');

micasplot
set(gcf, 'Position', [160, 1360, 800, 400]);
saveas(gcf, './Voutdiff_vs_DeltaZ.png');
