clear all
close all
clc

%% ==== SENSOR/READOUT PARAMS =====
N    = 10;
Cf   = 40e-15;      % feedback cap for readout: 40 fF
Catt = 0;           % attenuation cap
VTX  = 3;           % excitation amplitude (V)
VDD  = 3.3;
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

%% ==== Vout_diff (Vod) and SNR ====
% Vod = (VTX*(Cf+Catt)*(Csp-Csn)) / (Cf*(Csp+Csn+2*Cf+2*Catt))
Vod = (VTX * (Cf + Catt) .* (Csp - Csn)) ./ ...
      (Cf * (Csp + Csn + 2*Cf + 2*Catt));        % [V]

Vnoise_rms = 267e-6;                      % 267 uV RMS white noise

% Linear SNR (always positive)
SNR_lin = abs(Vod) ./ Vnoise_rms;         % |signal| / noise

% SNR in dB (can be negative when SNR_lin < 1)
SNR_dB  = 20*log10(SNR_lin);              % negative if signal < noise
SNR_dB(~isfinite(SNR_dB)) = NaN;          % clean up Inf/NaN

validSNR = SNR_dB(~isnan(SNR_dB));
minSNR   = min(validSNR);
maxSNR   = max(validSNR);

%% ==== Plot 1: Csp and Csn vs Delta-X for selected Delta-Z ==== 
figure('Color','w');
hold on;

num_steps = 5;
z_indices = round(linspace(1, length(deltaZ), num_steps));
colors = lines(num_steps);

hPlots  = [];
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

hSum   = gobjects(num_steps,1);
sumLeg = cell(1,num_steps);

for i = 1:num_steps
    idx  = z_indices(i);
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

%% ==== Plot 3: Vod across Delta-X and Delta-Z (2D map) + noise contour ====
figure('Color','w');

imagesc(deltaX*1e6, deltaZ*1e6, Vod);   % x: ΔX [µm], y: ΔZ [µm]
set(gca,'YDir','normal');               % so y axis goes from 0 up
xlabel('\DeltaX (µm)');
ylabel('\DeltaZ (µm)');
cb = colorbar;
ylabel(cb, 'V_{out\_diff} (V)');
title('V_{out\_diff} across \DeltaX and \DeltaZ');

hold on;
contour(deltaX*1e6, deltaZ*1e6, abs(Vod), ...
        [Vnoise_rms Vnoise_rms], 'k--', 'LineWidth', 1.5);

grid on;

micasplot
set(gcf, 'Position', [100, 100, 800, 500]);
saveas(gcf, './Vod_2D_DR.png');

%% ==== Plot 4: Vod vs Delta-X for selected Delta-Z slices + dashed-fill noise band + inset (below legend) ====
figure('Color','w');
hold on;

% ---- Make Plot 4 wider (room for legend + inset on the right)
set(gcf, 'Position', [120, 650, 1050, 520]);

% ---- Draw black dashed "fill" between the horizontal noise lines
x_um = deltaX * 1e6;
x_min = min(x_um);
x_max = max(x_um);

y1 =  Vnoise_rms;
y2 = -Vnoise_rms;

hBand = patch([x_min x_max x_max x_min], [y2 y2 y1 y1], ...
              'k', 'EdgeColor','none', 'FaceAlpha', 0.03);
set(hBand, 'HandleVisibility','off');

nFillLines = 22;
yFill = linspace(y2, y1, nFillLines);

nSeg   = 60;
dashOn = 0.60;
xEdges = linspace(x_min, x_max, nSeg+1);

for yy = yFill
    for k = 1:nSeg
        xa = xEdges(k);
        xb = xEdges(k+1);
        xb_on = xa + dashOn*(xb - xa);
        plot([xa xb_on], [yy yy], 'k-', 'LineWidth', 0.45, ...
             'HandleVisibility','off');
    end
end

% ---- Plot Vod curves
num_steps = 5;
z_indices = round(linspace(1, length(deltaZ), num_steps));
colors    = lines(num_steps);

hVod   = [];
legVod = {};

for i = 1:num_steps
    idx = z_indices(i);
    h = plot(x_um, Vod(idx,:), 'Color', colors(i,:), 'LineWidth', 1.6);
    hVod   = [hVod, h];
    legVod = [legVod, {sprintf('\\DeltaZ = %.1f µm', deltaZ(idx)*1e6)}];
end

yline(Vnoise_rms,  'k--', 'LineWidth', 1.5);
yline(-Vnoise_rms, 'k--', 'LineWidth', 1.5);

grid on
xlabel('{\Delta}X (µm)');
ylabel('V_{out\_diff} (V)');
legend(hVod, legVod, 'Location','bestoutside');

% -------- Inset: zoom where |Vod| < Vnoise_rms (place below legend area) --------
axMain = gca;

mid_idx  = z_indices(round(num_steps/2));
Vod_mid  = Vod(mid_idx,:);
mask_low = abs(Vod_mid) < Vnoise_rms;

if any(mask_low)
    x_zoom = deltaX(mask_low) * 1e6;
    y_zoom = Vod_mid(mask_low);

    x_minz = min(x_zoom);
    x_maxz = max(x_zoom);

    x_pad = 0.10 * (x_maxz - x_minz);
    if x_pad == 0, x_pad = 1; end

    y_minz = min(y_zoom);
    y_maxz = max(y_zoom);
    y_pad  = 0.25 * (y_maxz - y_minz + eps);

    mainPos  = get(axMain, 'Position');

    insetW   = 0.34 * mainPos(3);
    insetH   = 0.34 * mainPos(4);

    insetL   = mainPos(1) + 1.02*mainPos(3);
    insetB   = mainPos(2) + 0.10*mainPos(4);

    axInset  = axes('Position', [insetL insetB insetW insetH]);
    box(axInset, 'on');
    hold(axInset, 'on');

    for i = 1:num_steps
        idx = z_indices(i);
        plot(axInset, x_um, Vod(idx,:), 'Color', colors(i,:), 'LineWidth', 1.0);
    end

    yline(axInset, Vnoise_rms,  'k--', 'LineWidth', 1.0);
    yline(axInset, -Vnoise_rms, 'k--', 'LineWidth', 1.0);

    xlim(axInset, [x_minz - x_pad, x_maxz + x_pad]);
    ylim(axInset, [y_minz - y_pad, y_maxz + y_pad]);

    set(axInset, 'FontSize', 8);
    xlabel(axInset, '\DeltaX (µm)', 'FontSize', 8);
    ylabel(axInset, 'V_{out\_diff} (V)', 'FontSize', 8);

    % ---- Raise the title more (higher than previous)
    t = title(axInset, '|V_{od}| < V_{noise}', 'FontSize', 8);
    set(t, 'Units','normalized');
    tpos = get(t, 'Position');
    tpos(2) = tpos(2) + 0.18;   % <-- moved higher (was +0.08)
    set(t, 'Position', tpos);
end

micasplot
saveas(gcf, './Vod_vs_DX_DR.png');

%% ==== Plot 5: SNR (dB) across Delta-X and Delta-Z ====
figure('Color','w');

imagesc(deltaX*1e6, deltaZ*1e6, SNR_dB);
set(gca,'YDir','normal');
xlabel('\DeltaX (µm)');
ylabel('\DeltaZ (µm)');
cb = colorbar;
ylabel(cb, 'SNR (dB)');
title('SNR (dB) across \DeltaX and \DeltaZ, V_{noise} = 267 \muV_{RMS}');
grid on;

caxis([minSNR maxSNR]);

micasplot
set(gcf, 'Position', [150, 150, 800, 500]);
saveas(gcf, './SNR_2D_DR.png');

%% ==== Plot 6: SNR (dB) vs Delta-X for selected Delta-Z slices ====
figure('Color','w');
hold on;

hSNR   = [];
legSNR = {};

mask_pos = deltaX > 0;
dx_pos_um = deltaX(mask_pos) * 1e6;

for i = 1:num_steps
    idx = z_indices(i);
    h = plot(dx_pos_um, SNR_dB(idx, mask_pos), ...
        'Color', colors(i,:), 'LineWidth', 1.6);
    hSNR   = [hSNR, h];
    legSNR = [legSNR, {sprintf('\\DeltaZ = %.1f µm', deltaZ(idx)*1e6)}];
end

yline(0, 'k--', 'LineWidth', 1.5);

grid on
xlabel('{\Delta}X (µm)');
ylabel('SNR (dB)');
legend(hSNR, legSNR, 'Location','bestoutside');
set(gca,'XScale','log')

micasplot
set(gcf, 'Position', [180, 700, 800, 400]);
saveas(gcf, './SNR_vs_DX_dB_DR.png');


