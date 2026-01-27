%% Sweep sensor count and sensor pitch (plot only what exists) - SEPARATE SNN & SVM FIGURES
% - ONLY vertical Merkel reference lines (NO horizontal red dashed lines anywhere)
% - "Merkel cell pitch" label text is 4x larger
% - Vertical Merkel dashed line is 4x thicker
% - Legends for sensor pitch include "mm"
% - NO figi increment (every run overwrites same figure numbers and filenames)
% - Always runs micasplot before saving each figure
clearvars -except figi
close all

figi = 69;          % base figure number (kept constant)
kfolds = 6;

% ---------------- USER OPTIONS: FILTERING ----------------
mr_cnt_keep     = [];   % [] = keep all
spatialres_keep = [];   % [] = keep all
mr_cnt_exclude     = [];  % [] = exclude none
spatialres_exclude = [];  % [] = exclude none
% ---------------------------------------------------------

% Full grids
spatialres_array_full = [0.025, 0.05, 0.1, 0.2, 0.4, 1, 2, 5];   % sensor pitch (mm)
mr_cnt_array_full     = [3*3, 8*8, 14*14, 16*16, 18*18];         % sensor count

baseRoot = '../spikeclassifier_2/';

% ---------------- Apply filtering ----------------
mr_cnt_array = mr_cnt_array_full;
spatialres_array = spatialres_array_full;

if ~isempty(mr_cnt_keep)
    mr_cnt_array = mr_cnt_array(ismember(mr_cnt_array, mr_cnt_keep));
end
if ~isempty(mr_cnt_exclude)
    mr_cnt_array = mr_cnt_array(~ismember(mr_cnt_array, mr_cnt_exclude));
end

if ~isempty(spatialres_keep)
    spatialres_array = spatialres_array(ismember(spatialres_array, spatialres_keep));
end
if ~isempty(spatialres_exclude)
    spatialres_array = spatialres_array(~ismember(spatialres_array, spatialres_exclude));
end

if isempty(mr_cnt_array) || isempty(spatialres_array)
    error('After filtering, sensor count or sensor pitch array is empty. Adjust keep/exclude options.');
end

% Merkel cell pitch reference (for sensor pitch axis)
ref_sr = 0.4;
use_merkel_vertical = any(abs(spatialres_array - ref_sr) < 1e-12); % only draw line if in current x-grid

%% Storage (NaN means missing)
snn_mean = nan(numel(mr_cnt_array), numel(spatialres_array));
snn_std  = nan(numel(mr_cnt_array), numel(spatialres_array));

svm_mean = nan(numel(mr_cnt_array), numel(spatialres_array));
svm_std  = nan(numel(mr_cnt_array), numel(spatialres_array));

svm_bestX = nan(numel(mr_cnt_array), numel(spatialres_array));

% Store raw points per condition (optional use later)
snn_vals = cell(numel(mr_cnt_array), numel(spatialres_array)); % each cell: vector of fold accuracies (0..1)
svm_vals = cell(numel(mr_cnt_array), numel(spatialres_array)); % each cell: bestVec (0..1)

% --- helper to format sensor pitch exactly like folder names (e.g., 0.10 -> 0.1, 2.000 -> 2)
fmt_sr = @(x) regexprep(regexprep(sprintf('%.6f', x), '0+$', ''), '\.$', '');

%% ---------------- Collect results (skip missing folders/files) ----------------
for im = 1:numel(mr_cnt_array)
    mr_cnt = mr_cnt_array(im);

    % Python naming logic: x_size fixed 16, y_size depends on mr_cnt
    x_size = 16;
    y_size = ceil(mr_cnt / x_size);

    for is = 1:numel(spatialres_array)
        spatialres = spatialres_array(is);
        sr_str = fmt_sr(spatialres);

        %% ---------- SNN ----------
        snnFolder = sprintf([ ...
            'SNN_ep-150-isize-2-16-%d-hsize-256-bsize-32-lr-0.001-' ...
            'tmax_ms-400.0-Ts-10.0-k-6-dsetPath-ALL_N_lcadc-0-' ...
            'uV2-0-uV2-0-uV2-0-mult-0-uV2-3-multitrial-' ...
            'Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to40-1to40-' ...
            'mr_cnt-%d-spatialres-%s_sorted' ...
            ], y_size, mr_cnt, sr_str);

        snnPath = [baseRoot snnFolder filesep];
        snnFile = fullfile(snnPath, 'results_max.txt');

        if isfolder(snnPath) && isfile(snnFile)
            r_snn = parse_result(snnPath); % vector (0..1)
            if ~isempty(r_snn) && all(~isnan(r_snn))
                snn_vals{im,is} = r_snn;
                snn_mean(im,is) = mean(r_snn);
                snn_std(im,is)  = std(r_snn);
            end
        end

        %% ---------- SVM (best-0..5, choose best by mean) ----------
        bestMean = -inf;
        bestVec  = [];
        bestIdx  = NaN;

        for b = 0:(kfolds-1)
            svmFolder = sprintf([ ...
                'SVM_svmep-5-svm_wsize-40-sensor_sweep-lcadc-0-uV2-0-uV2-0-' ...
                'uV2-0-mult-0-uV2-3-mr_cnt-%d-spatialres-%s-best-%d' ...
                ], mr_cnt, sr_str, b);

            svmPath = [baseRoot svmFolder filesep];
            svmFile = fullfile(svmPath, 'results_svm_linear_max.txt');

            if isfolder(svmPath) && isfile(svmFile)
                v = parse_result_svm(svmPath); % vector (0..1)
                if ~isempty(v) && all(~isnan(v))
                    m = mean(v);
                    if m > bestMean
                        bestMean = m;
                        bestVec  = v;
                        bestIdx  = b;
                    end
                end
            end
        end

        if ~isempty(bestVec)
            svm_vals{im,is}  = bestVec;
            svm_bestX(im,is) = bestIdx;
            svm_mean(im,is)  = mean(bestVec);
            svm_std(im,is)   = std(bestVec);
        end
    end
end

%% ---------------- Normalize vs sensor pitch = 0.4 (Merkel cell pitch) ----------------
[~, ref_idx] = min(abs(spatialres_array - ref_sr));
can_normalize = (abs(spatialres_array(ref_idx) - ref_sr) < 1e-12);

snn_mean_norm = nan(size(snn_mean));
snn_std_norm  = nan(size(snn_std));

svm_mean_norm = nan(size(svm_mean));
svm_std_norm  = nan(size(svm_std));

if can_normalize
    for im = 1:numel(mr_cnt_array)
        snn_ref = snn_mean(im, ref_idx);
        svm_ref = svm_mean(im, ref_idx);

        if ~isnan(snn_ref) && snn_ref ~= 0
            snn_mean_norm(im,:) = snn_mean(im,:) ./ snn_ref;
            snn_std_norm(im,:)  = snn_std(im,:)  ./ snn_ref;
        end

        if ~isnan(svm_ref) && svm_ref ~= 0
            svm_mean_norm(im,:) = svm_mean(im,:) ./ svm_ref;
            svm_std_norm(im,:)  = svm_std(im,:)  ./ svm_ref;
        end
    end
else
    warning('sensor pitch=0.4 not included after filtering; normalized plots will be all-NaN and will say no results.');
end

%% ---------------- Plot helpers ----------------
% Vertical Merkel line with manual big label (compatible with older MATLAB).
add_merkel_vertical = @() add_merkel_vertical_manual(ref_sr);

% Always run micasplot + standard save sizing
width = 800; height = 700;
save_micas = @(fname_base) local_save_with_micasplot(fname_base, width, height);

% Multi-line ylabel strings (use newline, never "\n")
yl_raw_snn  = ['Classification Accuracy' newline '(SNN)'];
yl_raw_svm  = ['Classification Accuracy' newline '(SNN+SVM)'];
yl_norm_snn = ['Normalized accuracy to accuracy' newline ...
               'at Merkel cell pitch' newline '(SNN)'];
yl_norm_svm = ['Normalized accuracy to accuracy' newline ...
               'at Merkel cell pitch' newline '(SNN+SVM)'];

%% =======================================================================
%% 1) sensor pitch-sweep plots (with Merkel vertical line if applicable)
%% =======================================================================

% Figure 1: SNN raw - errorbars
figure(figi); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = snn_mean(im,:);
    e = snn_std(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_raw_snn);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SNN-raw-errorbar');

% Figure 2: SVM raw - errorbars
figure(figi+1); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = svm_mean(im,:);
    e = svm_std(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_raw_svm);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SVM-raw-errorbar');

% Figure 3: SNN norm - errorbars
figure(figi+2); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = snn_mean_norm(im,:);
    e = snn_std_norm(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_norm_snn);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SNN-norm-errorbar');

% Figure 4: SVM norm - errorbars
figure(figi+3); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = svm_mean_norm(im,:);
    e = svm_std_norm(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_norm_svm);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SVM-norm-errorbar');

%% =======================================================================
%% 2) Shaded variance envelope plots (color matches mean; legend fixed)
%% =======================================================================

% Figure 5: SNN raw - shaded variance
figure(figi+4); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = snn_mean(im,:);
    e = snn_std(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_raw_snn);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SNN-raw-band');

% Figure 6: SVM raw - shaded variance
figure(figi+5); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = svm_mean(im,:);
    e = svm_std(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_raw_svm);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SVM-raw-band');

% Figure 7: SNN norm - shaded variance
figure(figi+6); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = snn_mean_norm(im,:);
    e = snn_std_norm(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_norm_snn);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SNN-norm-band');

% Figure 8: SVM norm - shaded variance
figure(figi+7); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = svm_mean_norm(im,:);
    e = svm_std_norm(im,:);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor count=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end
set(gca,'XScale','log'); grid on;
xlabel('sensor pitch'); ylabel(yl_norm_svm);
if use_merkel_vertical, add_merkel_vertical(); end
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorPitch_sweep-SVM-norm-band');

%% =======================================================================
%% 3) sensor count sweep plots for each sensor pitch (NO Merkel vertical line) + bands
%%     Legend entries: "sensor pitch = X mm"
%% =======================================================================

% Figure 9: SNN raw - sweep sensor count (by sensor pitch)
figure(figi+8); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for is = 1:numel(spatialres_array)
    sr = spatialres_array(is);
    x  = mr_cnt_array(:);
    y  = snn_mean(:,is);
    e  = snn_std(:,is);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor pitch = %g mm', sr); %#ok<AGROW>
    end
end
grid on;
xlabel('sensor count'); ylabel(yl_raw_snn);
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorCount_sweep-SNN-raw-band');

% Figure 10: SVM raw - sweep sensor count (by sensor pitch)
figure(figi+9); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for is = 1:numel(spatialres_array)
    sr = spatialres_array(is);
    x  = mr_cnt_array(:);
    y  = svm_mean(:,is);
    e  = svm_std(:,is);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor pitch = %g mm', sr); %#ok<AGROW>
    end
end
grid on;
xlabel('sensor count'); ylabel(yl_raw_svm);
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorCount_sweep-SVM-raw-band');

% Figure 11: SNN norm - sweep sensor count (by sensor pitch)
figure(figi+10); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for is = 1:numel(spatialres_array)
    sr = spatialres_array(is);
    x  = mr_cnt_array(:);
    y  = snn_mean_norm(:,is);
    e  = snn_std_norm(:,is);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor pitch = %g mm', sr); %#ok<AGROW>
    end
end
grid on;
xlabel('sensor count'); ylabel(yl_norm_snn);
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SNN results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorCount_sweep-SNN-norm-band');

% Figure 12: SVM norm - sweep sensor count (by sensor pitch)
figure(figi+11); clf; hold on;
hLeg = gobjects(0); legTxt = {};
for is = 1:numel(spatialres_array)
    sr = spatialres_array(is);
    x  = mr_cnt_array(:);
    y  = svm_mean_norm(:,is);
    e  = svm_std_norm(:,is);
    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        h = plot(x(mask), y(mask), '-o', 'LineWidth', 2);
        if nnz(mask) >= 2
            band_from_line(x(mask), y(mask), e(mask), h);
            uistack(h,'top');
        end
        hLeg(end+1) = h; %#ok<AGROW>
        legTxt{end+1} = sprintf('sensor pitch = %g mm', sr); %#ok<AGROW>
    end
end
grid on;
xlabel('sensor count'); ylabel(yl_norm_svm);
if ~isempty(hLeg), legend(hLeg, legTxt, 'Location','best');
else, text(0.5,0.5,'No SVM results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center'); end
save_micas('./sensorCount_sweep-SVM-norm-band');

%% Optional: print chosen best-X table (NaN = missing)
disp('Chosen best-X (rows=sensor count, cols=sensor pitch) [NaN=missing]:');
disp(array2table(svm_bestX, ...
    'RowNames', cellfun(@(x) sprintf('sc_%s',x), string(mr_cnt_array), 'UniformOutput', false), ...
    'VariableNames', cellfun(@(x) matlab.lang.makeValidName(sprintf('sp_%s',x)), ...
                             string(spatialres_array), 'UniformOutput', false)));

%% ---------------- Local functions (MUST be at end of script) ----------------
function results_max = parse_result(resultsPath)
    results_max = [];

    fname = fullfile(resultsPath, 'results_max.txt');
    fid = fopen(fname, 'r');
    if fid < 0
        results_max = [];
        return;
    end

    tline = fgets(fid);
    while ischar(tline) && ~isempty(strtrim(tline))
        toks = split(tline, ':');
        if numel(toks) >= 2
            val = str2double(strtrim(toks{2}));
            results_max = [results_max; val]; %#ok<AGROW>
        end
        if numel(results_max) >= 6
            break
        end
        tline = fgets(fid);
    end

    fclose(fid);
    results_max = results_max ./ 100;
    results_max = results_max(~isnan(results_max));
end

function results_linear_max = parse_result_svm(resultsPath)
    results_linear_max = [];

    fname = fullfile(resultsPath, 'results_svm_linear_max.txt');
    fid = fopen(fname, 'r');
    if fid < 0
        results_linear_max = [];
        return;
    end

    tline = fgets(fid);
    while ischar(tline) && ~isempty(strtrim(tline))
        toks = split(tline, ':');
        if numel(toks) >= 2
            val = str2double(strtrim(toks{2}));
            results_linear_max = [results_linear_max; val]; %#ok<AGROW>
        end
        if numel(results_linear_max) >= 10
            break
        end
        tline = fgets(fid);
    end

    fclose(fid);
    results_linear_max = results_linear_max ./ 100;
    results_linear_max = results_linear_max(~isnan(results_linear_max));
end

function hb = band_from_line(x, y, e, hline)
    c = hline.Color;
    hb = fill([x(:); flipud(x(:))], [ (y(:)+e(:)); flipud(y(:)-e(:)) ], ...
              c, 'FaceAlpha',0.18, 'EdgeColor','none');
end

function add_merkel_vertical_manual(x0)
    % Draw the dashed vertical line (4x thicker than default)
    xl = xline(x0, 'r--', 'LineWidth', 4);
    xl.HandleVisibility = 'off';

    % Make a big label via text() positioned near the top of axes
    ax = gca;
    yl = ylim(ax);
    y_text = yl(2) - 0.06*(yl(2)-yl(1)); % slightly below top

    % On log-x axes, offset multiplicatively; on linear, offset additively
    if strcmpi(ax.XScale,'log')
        x_text = x0 * 1.08;
    else
        xlimv = xlim(ax);
        x_text = x0 + 0.02*(xlimv(2)-xlimv(1));
    end

    % Use axes font size as baseline, 4x larger
    fs = get(ax, 'FontSize');
    text(x_text, y_text, 'Merkel cell pitch', ...
        'Color','r', 'FontSize', fs*4, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Clipping','on');
end

function local_save_with_micasplot(fname_base, width, height)
    % Always run micasplot (assumed available on path)
    micasplot;

    % Standard size
    set(gcf, 'Position', [100, 100, width, height]);

    % Save both png + fig (overwrites each run)
    saveas(gcf, [fname_base '.png']);
    saveas(gcf, [fname_base '.fig']);
end


