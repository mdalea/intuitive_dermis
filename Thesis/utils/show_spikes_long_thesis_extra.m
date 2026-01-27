%% USED TO PLOT SPIKES DIRECTLY GATHERED THROUGH MATLAB (read_serial_spikes.m)
% ALWAYS SAVES ALL FIGURES/IMAGES
% close all

clearvars
clc

folderPath = pwd;
name='sample_6_adc_ibias_afe_5n_ibias_adapt_20n_islew_0p1u_vtrig_0p5_vin_6000mV_fin_1_3v_3v_3p5v-2025-08-15_02-59-27_sorted';
trial=1;
chunks=0;  % 0 if dataset without the -x

if chunks==0
    append = '';
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/',name,'/',num2str(trial),append,'.bs2'];
TD = Read_Ndataset_longTs(filename);

%% ---- PARAMETERS (requested) ----
fin      = 1;      % Hz
T        = 1/fin;  % s
t0       = 73.7;   % s (starting time)
vamp_vpp = 6;      % Vpp
vofs     = 0;      % V
maxPeriodsToPlot = 10;   % limit to first 10 periods only

% convert timestamps to seconds (TD.ts in us in your plotting code)
t = double(TD.ts)/1e6;
t_end = max(t);

% keep only *full* non-truncated periods from t0
nPeriodsAvail = floor((t_end - t0)/T);
if nPeriodsAvail <= 0
    error('No full periods available after t0=%.3f s (t_end=%.3f s).', t0, t_end);
end

% limit to first 10 periods
nPeriods  = min(nPeriodsAvail, maxPeriodsToPlot);
t_fit_end = t0 + nPeriods*T;

%% ---- SAVE DIRECTORY (always save figures/images) ----
outDir = fullfile(folderPath,'out_figs',name);
if ~exist(outDir,'dir'); mkdir(outDir); end

%% ---- Raster plot (all spikes) ----
fig_raster = figure('Name','Raster (all spikes)','Color','w');
hold on;

for spikeCount = 1:length(TD.ts)
    ylin = xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16);
    if TD.p(spikeCount)==1  % ON spike
        plot([t(spikeCount) t(spikeCount)], [ylin-0.4 ylin+0.4], 'k');
    else                   % OFF spike
        plot([t(spikeCount) t(spikeCount)], [ylin-0.4 ylin+0.4], 'Color','red');
    end
end

grid on
xlabel('Time (s)');
ylabel('Taxel Number');
title([folderPath,'/out/',name,'/',num2str(trial),append,'.bs2']);
micasplot

savefig(fig_raster, fullfile(outDir, sprintf('raster_all_trial%d%s.fig',trial,append)));
exportgraphics(fig_raster, fullfile(outDir, sprintf('raster_all_trial%d%s.png',trial,append)), 'Resolution', 300);

%% ---- Split ON/OFF ----
TD_on  = struct();
TD_off = struct();

TD_on.x  = TD.x(TD.p==1);
TD_on.y  = TD.y(TD.p==1);
TD_on.p  = TD.p(TD.p==1);
TD_on.ts = TD.ts(TD.p==1);

TD_off.x  = TD.x(TD.p==2);
TD_off.y  = TD.y(TD.p==2);
TD_off.p  = TD.p(TD.p==2);
TD_off.ts = TD.ts(TD.p==2);

t_on  = double(TD_on.ts)/1e6;
t_off = double(TD_off.ts)/1e6;

%% ---- Count ON/OFF spikes per (first 10) full periods and histogram ----
edges = (0:1:200); % adjust if you expect >200 spikes/period

onCounts  = zeros(nPeriods,1);
offCounts = zeros(nPeriods,1);

for k = 1:nPeriods
    a = t0 + (k-1)*T;
    b = t0 + k*T;
    onCounts(k)  = sum(t_on  >= a & t_on  < b);
    offCounts(k) = sum(t_off >= a & t_off < b);
end

totalOn  = sum(onCounts);
totalOff = sum(offCounts);

fprintf('Full periods USED (limited): %d (from %.3f s to %.3f s)\n', nPeriods, t0, t_fit_end);
fprintf('Total ON spikes in USED periods:  %d\n', totalOn);
fprintf('Total OFF spikes in USED periods: %d\n', totalOff);

fig_hist = figure('Name','Spike count histogram (per period, first 10)','Color','w');
hold on;

h1 = histogram(onCounts,  edges, 'DisplayStyle','stairs', 'LineWidth', 2, 'EdgeColor','k');
h2 = histogram(offCounts, edges, 'DisplayStyle','stairs', 'LineWidth', 2, 'EdgeColor','r');

grid on;
xlabel(sprintf('Spike count per %.3f s period (fin=%.3f Hz)', T, fin));
ylabel('Number of periods');
title(sprintf('First %d full periods only: [%.3f, %.3f] s', nPeriods, t0, t_fit_end));
legend([h1 h2], {'ON','OFF'}, 'Location','best');
micasplot

savefig(fig_hist, fullfile(outDir, sprintf('hist_counts_per_period_first%d_trial%d%s.fig',nPeriods,trial,append)));
exportgraphics(fig_hist, fullfile(outDir, sprintf('hist_counts_per_period_first%d_trial%d%s.png',nPeriods,trial,append)), 'Resolution', 300);

%% ---- LSB/2 metrics PER PERIOD (as requested) ----
% For each period k:
%   LSB/2(k)  = Vpp / N_on(k)
%   -LSB/2(k) = Vpp / N_off(k)
% Then report averages over the plotted periods.
LSB_over2_perPeriod_pos = nan(nPeriods,1);
LSB_over2_perPeriod_neg = nan(nPeriods,1);

for k = 1:nPeriods
    if onCounts(k)  > 0, LSB_over2_perPeriod_pos(k) = vamp_vpp / onCounts(k);  end
    if offCounts(k) > 0, LSB_over2_perPeriod_neg(k) = vamp_vpp / offCounts(k); end
end

avg_LSB_over2_pos = mean(LSB_over2_perPeriod_pos(~isnan(LSB_over2_perPeriod_pos)));
avg_LSB_over2_neg = mean(LSB_over2_perPeriod_neg(~isnan(LSB_over2_perPeriod_neg)));

fprintf('avg (LSB/2)  over first %d periods = mean(Vpp./N_on(k))  = %.6g V\n', nPeriods, avg_LSB_over2_pos);
fprintf('avg (-LSB/2) over first %d periods = mean(Vpp./N_off(k)) = %.6g V\n', nPeriods, avg_LSB_over2_neg);

%% ---- Figure: BIG sinusoid + SMALL aligned ON/OFF spikes (FIRST 10 PERIODS ONLY) ----
% restrict spikes to [t0, t_fit_end)
mask_on  = (t_on  >= t0) & (t_on  < t_fit_end);
mask_off = (t_off >= t0) & (t_off < t_fit_end);

t_on_win  = t_on(mask_on);
t_off_win = t_off(mask_off);

% waveform timebase for plotting
fs_plot = 5000; % Hz (just for smooth display)
t_wave = linspace(t0, t_fit_end, max(2, round((t_fit_end - t0)*fs_plot)));

vpk = vamp_vpp/2;
v   = vpk * sin(2*pi*fin*t_wave) + vofs;

fig_overlay = figure('Name','Sinusoid (big) + aligned ON/OFF spikes (small) - first 10 periods','Color','w');

%subplot(1,2,1);
hold on;

% BIG sinusoid
yyaxis left
plot(t_wave, v, 'LineWidth', 4);
ylabel('Vin (V)');
ylim([vofs - 1.2*vpk, vofs + 1.2*vpk]);
grid on;

% SMALL spikes
yyaxis right
spkAmp = 0.15;
spkLW  = 0.75;

for i = 1:numel(t_on_win)
    plot([t_on_win(i) t_on_win(i)], [0 spkAmp], 'k', 'LineWidth', spkLW);
end
for i = 1:numel(t_off_win)
    plot([t_off_win(i) t_off_win(i)], [0 -spkAmp], 'Color','red', 'LineWidth', spkLW);
end

ylim([-0.35 0.35]);
yticks([-spkAmp 0 spkAmp]);
ylabel('Spikes (polarity)');
xlim([t0 t_fit_end]);
xlabel('Time (s)');
%title(sprintf('fin=%.3f Hz, Vpp=%g V, first %d periods', fin, vamp_vpp, nPeriods));
%micasplot

% Right subplot empty (no text), as requested
%subplot(1,2,2);
axis off;

savefig(fig_overlay, fullfile(outDir, sprintf('sinusoid_spikes_overlay_first%d_trial%d%s.fig',nPeriods,trial,append)));
exportgraphics(fig_overlay, fullfile(outDir, sprintf('sinusoid_spikes_overlay_first%d_trial%d%s.png',nPeriods,trial,append)), 'Resolution', 300);
