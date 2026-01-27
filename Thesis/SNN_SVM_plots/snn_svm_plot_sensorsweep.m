%% Sweep mr_cnt and spatialres (plot only what exists) - SEPARATE SNN & SVM FIGURES
% Normalize accuracy versus accuracy at spatialres = 0.4 (Merkel cell pitch)
clearvars -except figi
close all

figi = 69;
kfolds = 6;

spatialres_array = [0.025, 0.05, 0.1, 0.2, 0.4, 1, 2, 5];
mr_cnt_array     = [3*3, 8*8, 14*14, 16*16, 18*18];

baseRoot = '../spikeclassifier/';

% Storage (NaN means missing)
snn_mean = nan(numel(mr_cnt_array), numel(spatialres_array));
snn_std  = nan(numel(mr_cnt_array), numel(spatialres_array));

svm_mean = nan(numel(mr_cnt_array), numel(spatialres_array));
svm_std  = nan(numel(mr_cnt_array), numel(spatialres_array));

svm_bestX = nan(numel(mr_cnt_array), numel(spatialres_array));

% --- helper to format spatialres exactly like folder names (e.g., 0.10 -> 0.1, 2.000 -> 2)
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
            r_snn = parse_result(snnPath); % returns vector
            if ~isempty(r_snn) && all(~isnan(r_snn))
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
                v = parse_result_svm(svmPath);
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
            svm_bestX(im,is) = bestIdx;
            svm_mean(im,is)  = mean(bestVec);
            svm_std(im,is)   = std(bestVec);
        end
    end
end

%% ---------------- Normalize vs spatialres = 0.4 (Merkel cell pitch) ----------------
ref_sr = 0.4;
[~, ref_idx] = min(abs(spatialres_array - ref_sr)); % robust to float representation

% normalized arrays (same shape)
snn_mean_norm = nan(size(snn_mean));
snn_std_norm  = nan(size(snn_std));

svm_mean_norm = nan(size(svm_mean));
svm_std_norm  = nan(size(svm_std));

for im = 1:numel(mr_cnt_array)
    snn_ref = snn_mean(im, ref_idx);
    svm_ref = svm_mean(im, ref_idx);

    % normalize mean and std by the reference mean (ratio normalization)
    if ~isnan(snn_ref) && snn_ref ~= 0
        snn_mean_norm(im,:) = snn_mean(im,:) ./ snn_ref;
        snn_std_norm(im,:)  = snn_std(im,:)  ./ snn_ref;
    end

    if ~isnan(svm_ref) && svm_ref ~= 0
        svm_mean_norm(im,:) = svm_mean(im,:) ./ svm_ref;
        svm_std_norm(im,:)  = svm_std(im,:)  ./ svm_ref;
    end
end

%% ---------------- Plot ONLY existing points (SEPARATE FIGURES) ----------------

% =========================
% Figure: SNN (normalized)
% =========================
figure(figi); clf; hold on;

legSNN = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = snn_mean_norm(im,:);
    e = snn_std_norm(im,:);

    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        legSNN{end+1} = sprintf('mr\\_cnt=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end

set(gca,'XScale','log');
grid on;
xlabel('spatialres');
ylabel('Normalized accuracy to accuracy at Merkel cell pitch');

if ~isempty(legSNN)
    legend(legSNN, 'Location','best');
else
    text(0.5,0.5,'No SNN results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center');
end

set(gcf,'Position',[100 100 650 450]);
saveas(gcf, sprintf('mrCnt_spatialRes_sweep-SNN-normMerkelPitch%.3g-%d.png', spatialres_array(ref_idx), figi));
saveas(gcf, sprintf('mrCnt_spatialRes_sweep-SNN-normMerkelPitch%.3g-%d.fig', spatialres_array(ref_idx), figi));

% =========================
% Figure: SVM (normalized)
% =========================
figure(figi+1); clf; hold on;

legSVM = {};
for im = 1:numel(mr_cnt_array)
    x = spatialres_array;
    y = svm_mean_norm(im,:);
    e = svm_std_norm(im,:);

    mask = ~isnan(y) & ~isnan(e);
    if nnz(mask) >= 1
        errorbar(x(mask), y(mask), e(mask), '-o', 'LineWidth', 2);
        legSVM{end+1} = sprintf('mr\\_cnt=%d', mr_cnt_array(im)); %#ok<AGROW>
    end
end

set(gca,'XScale','log');
grid on;
xlabel('spatialres');
ylabel('Normalized accuracy to accuracy at Merkel cell pitch');

if ~isempty(legSVM)
    legend(legSVM, 'Location','best');
else
    text(0.5,0.5,'No SVM results found (or missing 0.4 ref)','Units','normalized','HorizontalAlignment','center');
end

set(gcf,'Position',[800 100 650 450]);
saveas(gcf, sprintf('mrCnt_spatialRes_sweep-SVM-normMerkelPitch%.3g-%d.png', spatialres_array(ref_idx), figi+1));
saveas(gcf, sprintf('mrCnt_spatialRes_sweep-SVM-normMerkelPitch%.3g-%d.fig', spatialres_array(ref_idx), figi+1));

%% Optional: print chosen best-X table (NaN = missing)
disp('Chosen best-X (rows=mr_cnt, cols=spatialres) [NaN=missing]:');
disp(array2table(svm_bestX, ...
    'RowNames', cellfun(@(x) sprintf('mr_%s',x), string(mr_cnt_array), 'UniformOutput', false), ...
    'VariableNames', cellfun(@(x) matlab.lang.makeValidName(sprintf('sr_%s',x)), ...
                             string(spatialres_array), 'UniformOutput', false)));

% advance figure counter by 2 since we used two figures
figi = figi + 2;

%% ---------------- Local functions ----------------
function results_max = parse_result(resultsPath)
    results_max = [];

    fname = fullfile(resultsPath, 'results_max.txt');
    fid = fopen(fname, 'r');
    if fid < 0
        results_max = [];
        return;
    end

    % first line + next 5 lines = 6 vals (but tolerate shorter)
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

    % drop NaNs if any parsing failed
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

    % up to 10 vals (tolerate shorter)
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

    % drop NaNs if any parsing failed
    results_linear_max = results_linear_max(~isnan(results_linear_max));
end

