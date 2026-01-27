%% Sweep vs X (CRF size) @ fixed sensor count & pitch (SNN+SVM RESULTS ONLY)
clearvars -except figi
close all

figi = 69;

%% ---------------- USER OPTIONS ----------------
baseRoot = '../spikeclassifier_2/'; % where your sweep outputs live

mr_cnt_fixed = 14*14;   % 196
spatialres_fixed = 0.2;

X_array = [4 8 16 32 64 80 128 150];
Y_array = {'global','local'};
Z_array = {'unidist-max8','skewdist'};
K_array = 0:5;

outPrefix = './crf_sweep';

width = 900; height = 550;

X_bar = [32 128];
%% ---------------------------------------------

sp_str = strip_trailing_zeros(spatialres_fixed);

nX = numel(X_array);
nY = numel(Y_array);
nZ = numel(Z_array);

comboColors = lines(nY*nZ);
comboLabel = @(iy,iz) sprintf('%s, %s', Y_array{iy}, Z_array{iz});
comboColor = @(iy,iz) comboColors((iy-1)*nZ + iz, :);

%% ---------------- Storage ----------------
snn_mean = nan(nX,nY,nZ);
snn_std  = nan(nX,nY,nZ);
snn_vals = cell(nX,nY,nZ);

svm_mean  = nan(nX,nY,nZ);
svm_std   = nan(nX,nY,nZ);
svm_bestK = nan(nX,nY,nZ);
svm_vals  = cell(nX,nY,nZ);

%% ---------------- Collect results ----------------
for ix = 1:nX
    X = X_array(ix);

    for iy = 1:nY
        Y = Y_array{iy};

        for iz = 1:nZ
            Z = Z_array{iz};

            snnFolder = sprintf([ ...
                'SNN_dsetPath-ALL_N_lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-3-' ...
                'crf-%d-%s_crf-%s-' ...
                'multitrial-Kylberg_filt_6_scan_oversampled20x_actualdimscale--1to10-1to10-' ...
                'mr_cnt-%d-spatialres-%s_sorted' ...
                ], X, Y, Z, mr_cnt_fixed, sp_str);

            snnPath = fullfile(baseRoot, snnFolder);
            snnFile = fullfile(snnPath, 'results_max.txt');

            if isfolder(snnPath) && isfile(snnFile)
                v = parse_result_file(snnFile, 6);
                if ~isempty(v) && all(~isnan(v))
                    snn_vals{ix,iy,iz} = v;
                    snn_mean(ix,iy,iz) = mean(v);
                    snn_std(ix,iy,iz)  = std(v);
                end
            end

            bestMean = -inf;
            bestVec  = [];
            bestK    = NaN;

            for K = K_array
                svmFolder = sprintf([ ...
                    'SVM_svmep-5-svm_wsize-40-sensor_sweep_crf-' ...
                    'lcadc-0-uV2-0-uV2-0-uV2-0-mult-0-uV2-3-' ...
                    'crf-%d-%s_crf-%s-mr_cnt-%d-spatialres-%s-best-%d' ...
                    ], X, Y, Z, mr_cnt_fixed, sp_str, K);

                svmPath = fullfile(baseRoot, svmFolder);
                svmFile = fullfile(svmPath, 'results_svm_linear_max.txt');

                if isfolder(svmPath) && isfile(svmFile)
                    v = parse_result_file(svmFile, 10);
                    if ~isempty(v) && all(~isnan(v))
                        m = mean(v);
                        if m > bestMean
                            bestMean = m;
                            bestVec  = v;
                            bestK    = K;
                        end
                    end
                end
            end

            if ~isempty(bestVec)
                svm_vals{ix,iy,iz}  = bestVec;
                svm_bestK(ix,iy,iz) = bestK;
                svm_mean(ix,iy,iz)  = mean(bestVec);
                svm_std(ix,iy,iz)   = std(bestVec);
            end
        end
    end
end

%% ---------------- Labels ----------------
yl_snn = ['Classification Accuracy' newline '(SNN)'];
yl_svm = ['Classification Accuracy' newline '(SNN+SVM)'];

xlab = sprintf('Heminode count (raw taxel count = %d)', mr_cnt_fixed);

%% =======================================================================
%% 1) SNN vs X (errorbar)
%% =======================================================================
figure(figi); clf; hold on;
hLeg = gobjects(0); legTxt = {};

for iy = 1:nY
    for iz = 1:nZ
        x = X_array(:);
        y = squeeze(snn_mean(:,iy,iz));
        e = squeeze(snn_std(:,iy,iz));
        mask = ~isnan(y) & ~isnan(e);

        if nnz(mask) >= 1
            c = comboColor(iy,iz);
            h = errorbar(x(mask), y(mask), e(mask), '-o', ...
                'LineWidth', 2, 'Color', c, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            hLeg(end+1) = h; %#ok<AGROW>
            legTxt{end+1} = comboLabel(iy,iz); %#ok<AGROW>
        end
    end
end

grid on;
xlabel(xlab);
ylabel(yl_snn);
title(sprintf('SNN @ sensor pitch=%g mm', spatialres_fixed));

if ~isempty(hLeg)
    legend(hLeg, legTxt, 'Location','best');
else
    text(0.5,0.5,'No SNN results found','Units','normalized','HorizontalAlignment','center');
end

local_save_with_micasplot([outPrefix '_sweepX_SNN_fixed196_pitch0p2-errorbar'], width, height);

%% =======================================================================
%% 2) SNN vs X (shaded band)
%% =======================================================================
figure(figi+1); clf; hold on;
hLeg = gobjects(0); legTxt = {};

for iy = 1:nY
    for iz = 1:nZ
        x = X_array(:);
        y = squeeze(snn_mean(:,iy,iz));
        e = squeeze(snn_std(:,iy,iz));
        mask = ~isnan(y) & ~isnan(e);

        if nnz(mask) >= 1
            c = comboColor(iy,iz);
            h = plot(x(mask), y(mask), '-o', ...
                'LineWidth', 2, 'Color', c, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            if nnz(mask) >= 2
                band_from_line(x(mask), y(mask), e(mask), h);
                uistack(h,'top');
            end
            hLeg(end+1) = h; %#ok<AGROW>
            legTxt{end+1} = comboLabel(iy,iz); %#ok<AGROW>
        end
    end
end

grid on;
xlabel(xlab);
ylabel(yl_snn);
title(sprintf('SNN @ sensor pitch=%g mm', spatialres_fixed));

if ~isempty(hLeg)
    legend(hLeg, legTxt, 'Location','best');
else
    text(0.5,0.5,'No SNN results found','Units','normalized','HorizontalAlignment','center');
end

local_save_with_micasplot([outPrefix '_sweepX_SNN_fixed196_pitch0p2-band'], width, height);

%% =======================================================================
%% 3) SVM (SNN+SVM) vs X (errorbar, best-K)
%% =======================================================================
figure(figi+2); clf; hold on;
hLeg = gobjects(0); legTxt = {};

for iy = 1:nY
    for iz = 1:nZ
        x = X_array(:);
        y = squeeze(svm_mean(:,iy,iz));
        e = squeeze(svm_std(:,iy,iz));
        mask = ~isnan(y) & ~isnan(e);

        if nnz(mask) >= 1
            c = comboColor(iy,iz);
            h = errorbar(x(mask), y(mask), e(mask), '-o', ...
                'LineWidth', 2, 'Color', c, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            hLeg(end+1) = h; %#ok<AGROW>
            legTxt{end+1} = comboLabel(iy,iz); %#ok<AGROW>
        end
    end
end

grid on;
xlabel(xlab);
ylabel(yl_svm);
title(sprintf('SNN+SVM @ sensor pitch=%g mm (best K)', spatialres_fixed));

if ~isempty(hLeg)
    legend(hLeg, legTxt, 'Location','best');
else
    text(0.5,0.5,'No SVM results found','Units','normalized','HorizontalAlignment','center');
end

local_save_with_micasplot([outPrefix '_sweepX_SNNplusSVM_fixed196_pitch0p2-errorbar'], width, height);

%% =======================================================================
%% 4) SVM (SNN+SVM) vs X (shaded band, best-K)
%% =======================================================================
figure(figi+3); clf; hold on;
hLeg = gobjects(0); legTxt = {};

for iy = 1:nY
    for iz = 1:nZ
        x = X_array(:);
        y = squeeze(svm_mean(:,iy,iz));
        e = squeeze(svm_std(:,iy,iz));
        mask = ~isnan(y) & ~isnan(e);

        if nnz(mask) >= 1
            c = comboColor(iy,iz);
            h = plot(x(mask), y(mask), '-o', ...
                'LineWidth', 2, 'Color', c, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            if nnz(mask) >= 2
                band_from_line(x(mask), y(mask), e(mask), h);
                uistack(h,'top');
            end
            hLeg(end+1) = h; %#ok<AGROW>
            legTxt{end+1} = comboLabel(iy,iz); %#ok<AGROW>
        end
    end
end

grid on;
xlabel(xlab);
ylabel(yl_svm);
title(sprintf('SNN+SVM @ sensor pitch=%g mm (best K)', spatialres_fixed));

if ~isempty(hLeg)
    legend(hLeg, legTxt, 'Location','best');
else
    text(0.5,0.5,'No SVM results found','Units','normalized','HorizontalAlignment','center');
end

local_save_with_micasplot([outPrefix '_sweepX_SNNplusSVM_fixed196_pitch0p2-band'], width, height);

%% =======================================================================
%% 5) BAR PLOTS @ X in X_bar (same figure; colored to match vs-X)
%% =======================================================================
ix_bar = nan(size(X_bar));
for i = 1:numel(X_bar)
    ix_bar(i) = find(X_array == X_bar(i), 1);
    if isempty(ix_bar(i))
        error('X_array does not contain X=%d. Add it to X_array.', X_bar(i));
    end
end

iy_global = find(strcmp(Y_array,'global'), 1);
iy_local  = find(strcmp(Y_array,'local'),  1);
iz_uni    = find(strcmp(Z_array,'unidist-max8'), 1);
iz_skew   = find(strcmp(Z_array,'skewdist'),     1);

comboList = [ ...
    iy_global, iz_uni; ...
    iy_global, iz_skew; ...
    iy_local,  iz_uni; ...
    iy_local,  iz_skew ...
];

comboNames = { ...
    'global, unidist-max8', ...
    'global, skewdist', ...
    'local, unidist-max8', ...
    'local, skewdist' ...
};

nG = numel(X_bar);
nS = size(comboList,1);

snn_bar_mean = nan(nG, nS);
snn_bar_std  = nan(nG, nS);

svm_bar_mean = nan(nG, nS);
svm_bar_std  = nan(nG, nS);

for g = 1:nG
    ix = ix_bar(g);
    for s = 1:nS
        iy = comboList(s,1);
        iz = comboList(s,2);

        snn_bar_mean(g,s) = snn_mean(ix,iy,iz);
        snn_bar_std(g,s)  = snn_std(ix,iy,iz);

        svm_bar_mean(g,s) = svm_mean(ix,iy,iz);
        svm_bar_std(g,s)  = svm_std(ix,iy,iz);
    end
end

seriesColors = nan(nS,3);
for s = 1:nS
    seriesColors(s,:) = comboColor(comboList(s,1), comboList(s,2));
end

figure(figi+4); clf;
xcat = categorical(string(X_bar));
b = bar(xcat, snn_bar_mean, 'grouped'); hold on;
for s = 1:numel(b), b(s).FaceColor = seriesColors(s,:); end
for s = 1:numel(b)
    if isprop(b(s),'XEndPoints'), xx = b(s).XEndPoints; else, xx = 1:nG; end
    errorbar(xx, snn_bar_mean(:,s), snn_bar_std(:,s), 'k', 'LineStyle','none', 'LineWidth', 1.2);
end
grid on;
xlabel(xlab);
ylabel(yl_snn);
title('');
legend(comboNames, 'Location','best');
local_save_with_micasplot([outPrefix sprintf('_bar_SNN_X%s_mr%d_pitch%s', strjoin(string(X_bar),'_'), mr_cnt_fixed, sp_str)], width, height);

figure(figi+5); clf;
xcat = categorical(string(X_bar));
b = bar(xcat, svm_bar_mean, 'grouped'); hold on;
for s = 1:numel(b), b(s).FaceColor = seriesColors(s,:); end
for s = 1:numel(b)
    if isprop(b(s),'XEndPoints'), xx = b(s).XEndPoints; else, xx = 1:nG; end
    errorbar(xx, svm_bar_mean(:,s), svm_bar_std(:,s), 'k', 'LineStyle','none', 'LineWidth', 1.2);
end
grid on;
xlabel(xlab);
ylabel(yl_svm);
title('');
legend(comboNames, 'Location','best');
local_save_with_micasplot([outPrefix sprintf('_bar_SNNplusSVM_X%s_mr%d_pitch%s', strjoin(string(X_bar),'_'), mr_cnt_fixed, sp_str)], width, height);

%% ---------------- Local functions (MUST be at end of script) ----------------
function s = strip_trailing_zeros(x)
    s = sprintf('%.12f', x);
    s = regexprep(s, '0+$', '');
    s = regexprep(s, '\.$', '');
end

function vals01 = parse_result_file(fname, maxN)
    vals01 = [];
    fid = fopen(fname, 'r');
    if fid < 0, return; end
    tline = fgets(fid);
    while ischar(tline) && ~isempty(strtrim(tline))
        toks = split(tline, ':');
        if numel(toks) >= 2
            v = str2double(strtrim(toks{2}));
            vals01 = [vals01; v]; %#ok<AGROW>
        end
        if numel(vals01) >= maxN, break; end
        tline = fgets(fid);
    end
    fclose(fid);
    vals01 = vals01 ./ 100;
    vals01 = vals01(~isnan(vals01));
end

function hb = band_from_line(x, y, e, hline)
    c = hline.Color;
    hb = fill([x(:); flipud(x(:))], [ (y(:)+e(:)); flipud(y(:)-e(:)) ], ...
              c, 'FaceAlpha',0.18, 'EdgeColor','none');
end

function local_save_with_micasplot(fname_base, width, height)
    micasplot;
    set(gcf, 'Position', [100, 100, width, height]);
    saveas(gcf, [fname_base '.png']);
    saveas(gcf, [fname_base '.fig']);
end

