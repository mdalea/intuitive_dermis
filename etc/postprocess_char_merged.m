%% postprocess_merged.m
% MERGED script (ALL FOUR FAMILIES) + NO "family" prompt in sizing:
%
% Families processed:
%   A) flavors_p   : svt_p / lvt_p / hvt_p
%      - L CSV order preserved: [60nm, 500nm, 1u, 2u, 4u, 5u, 10u, 20u, 40u]
%      - Legend order preserved: same
%
%   B) devices_p   : 2p5v_p / 1p8v_p / 3p3v_p
%      - L CSV order preserved: [500nm, 1u, 2u, 4u, 5u, 10u, 20u, 40u]
%      - Legend order preserved: same
%
%   C) flavors     : svt / lvt / hvt
%      - L CSV order preserved: [500nm, 1u, 2u, 4u, 5u, 10u, 20u, 40u, 60nm]
%      - Legend order preserved: [60nm, 500nm, 1u, 2u, 4u, 5u, 10u, 20u, 40u]
%
%   D) devices     : 2p5v / 1p8v / 3p3v
%      - L CSV order preserved: [1u, 2u, 4u, 5u, 10u, 20u, 40u, 500nm]
%      - Legend order preserved: [500nm, 1u, 2u, 4u, 5u, 10u, 20u, 40u]
%
% COMMON requirements:
%   - normids uses ABSOLUTE VALUE everywhere (CSV output, plots, sizing)
%   - gmoverid vs normids uses LOG X-AXIS
%   - figures are kept OPEN
%   - robust plotting handle capture
%
% Sizing change per your request:
%   - user DOES NOT choose "family"
%   - user directly enters key: e.g. svt, svt_p, 2p5v, 2p5v_p, ...
%   - script auto-detects which family map contains that key

close all;
clear; clc;

%% --- User settings ---
inputDir  = pwd;                 % folder containing the input CSVs
outputDir = fullfile(pwd,'out'); % output folder
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

%% =========================
%  FAMILY DEFINITIONS
% =========================

% A) flavors_p
keys_A = {'svt_p','lvt_p','hvt_p'};
L_labels_csv_A    = {'60nm','500nm','1u','2u','4u','5u','10u','20u','40u'};
L_values_csv_A    = [60e-9, 500e-9, 1e-6, 2e-6, 4e-6, 5e-6, 10e-6, 20e-6, 40e-6];
L_labels_legend_A = {'60nm','500nm','1u','2u','4u','5u','10u','20u','40u'};
prefix_A = 'FLAVORP';

% B) devices_p
keys_B = {'2p5v_p','1p8v_p','3p3v_p'};
L_labels_csv_B    = {'500nm','1u','2u','4u','5u','10u','20u','40u'};
L_values_csv_B    = [500e-9, 1e-6, 2e-6, 4e-6, 5e-6, 10e-6, 20e-6, 40e-6];
L_labels_legend_B = {'500nm','1u','2u','4u','5u','10u','20u','40u'};
prefix_B = 'DEVICEP';

% C) flavors (no _p)
keys_C = {'svt','lvt','hvt'};
L_labels_csv_C    = {'500nm','1u','2u','4u','5u','10u','20u','40u','60nm'}; % CSV order preserved
L_values_csv_C    = [500e-9, 1e-6, 2e-6, 4e-6, 5e-6, 10e-6, 20e-6, 40e-6, 60e-9];
L_labels_legend_C = {'60nm','500nm','1u','2u','4u','5u','10u','20u','40u'}; % legend order preserved
prefix_C = 'FLAVOR';

% D) devices (NO _p)
keys_D = {'2p5v','1p8v','3p3v'};
L_labels_csv_D    = {'1u','2u','4u','5u','10u','20u','40u','500nm'}; % CSV order preserved
L_values_csv_D    = [1e-6, 2e-6, 4e-6, 5e-6, 10e-6, 20e-6, 40e-6, 500e-9];
L_labels_legend_D = {'500nm','1u','2u','4u','5u','10u','20u','40u'}; % legend order preserved
prefix_D = 'DEVICE';

%% Storage for sizing (Maps so "2p5v" keys are valid)
resultsA = containers.Map('KeyType','char','ValueType','any'); % flavor_p
resultsB = containers.Map('KeyType','char','ValueType','any'); % device_p
resultsC = containers.Map('KeyType','char','ValueType','any'); % flavor
resultsD = containers.Map('KeyType','char','ValueType','any'); % device

%% =========================
%  PROCESS ALL FAMILIES
% =========================
for i = 1:numel(keys_A)
    k = keys_A{i};
    resultsA(k) = processOneSet(inputDir, outputDir, k, prefix_A, L_labels_csv_A, L_values_csv_A, L_labels_legend_A);
end

for i = 1:numel(keys_B)
    k = keys_B{i};
    resultsB(k) = processOneSet(inputDir, outputDir, k, prefix_B, L_labels_csv_B, L_values_csv_B, L_labels_legend_B);
end

for i = 1:numel(keys_C)
    k = keys_C{i};
    resultsC(k) = processOneSet(inputDir, outputDir, k, prefix_C, L_labels_csv_C, L_values_csv_C, L_labels_legend_C);
end

for i = 1:numel(keys_D)
    k = keys_D{i};
    resultsD(k) = processOneSet(inputDir, outputDir, k, prefix_D, L_labels_csv_D, L_values_csv_D, L_labels_legend_D);
end

disp('All done. Outputs are in:');
disp(outputDir);

%% =========================
%  Interactive sizing helper (NO family selection)
% =========================
fprintf('\n=========================\n');
fprintf(' SIZING (gm/Id-based)\n');
fprintf('=========================\n');
fprintf('Type the device/flavor key directly (e.g. svt, svt_p, 2p5v, 2p5v_p).\n');
fprintf('L is FIXED once chosen; returns W_need for that L.\n');
fprintf('Uses abs(normids) everywhere.\n\n');

% Build a unified key list for the user
allKeys = [keys(resultsA), keys(resultsB), keys(resultsC), keys(resultsD)];
allKeys = unique(allKeys); % unique set
fprintf('Available keys:\n');
for k = 1:numel(allKeys)
    fprintf('  %s\n', allKeys{k});
end
fprintf('\n');

keepGoing = true;
while keepGoing
    try
        key = lower(strtrim(input('Enter device/flavor key: ', 's')));

        % Auto-detect which map contains the key (no family prompt)
        if isKey(resultsA, key)
            R = resultsA(key); famName = 'flavor_p';
        elseif isKey(resultsB, key)
            R = resultsB(key); famName = 'device_p';
        elseif isKey(resultsC, key)
            R = resultsC(key); famName = 'flavor';
        elseif isKey(resultsD, key)
            R = resultsD(key); famName = 'device';
        else
            error('Key "%s" not found. Please choose from the printed list.', key);
        end

        fprintf('\nSelected: %s (%s)\n', upper(key), upper(famName));

        fprintf('\nAvailable L options (legend order):\n');
        for k = 1:numel(R.L_labels)
            fprintf('  %s\n', R.L_labels{k});
        end

        L_sel = strtrim(input('Choose L label: ', 's'));
        idxL = find(strcmpi(R.L_labels, L_sel), 1, 'first');
        if isempty(idxL)
            error('L label "%s" not recognized.', L_sel);
        end
        L_chosen = R.L_m(idxL);

        gm_target = input('Enter desired gmoverid (gm/Id): ');
        if ~isscalar(gm_target) || ~isfinite(gm_target)
            error('gmoverid must be a finite scalar.');
        end

        gm_vec = R.gmoverid{idxL};
        ni_vec = abs(R.normids{idxL}); % safety
        normids_at_gm = interpNormIdsAtGm(gm_vec, ni_vec, gm_target);

        if ~isfinite(normids_at_gm) || normids_at_gm <= 0
            error(['Could not interpolate a positive |normalized Ids| at gmoverid=%.6g for L=%s. ' ...
                   'Check gm_target within curve range and data validity.'], gm_target, L_sel);
        end

        fprintf('\nInterpolated at L=%s:\n', L_sel);
        fprintf('  |normalized Ids| = %.6g (A / (W/L)) at gm/Id = %.6g\n', normids_at_gm, gm_target);

        Ids_need = input('Enter desired Ids_need magnitude (A, positive): ');
        if ~isscalar(Ids_need) || ~isfinite(Ids_need) || Ids_need <= 0
            error('Ids_need must be positive, finite scalar in Amps.');
        end

        WL_need = Ids_need / normids_at_gm;
        W_need  = WL_need * L_chosen;

        fprintf('\nSizing result:\n');
        fprintf('  Key         = %s\n', upper(key));
        fprintf('  L_need      = %s (%.6g m)\n', L_sel, L_chosen);
        fprintf('  (W/L)_need  = %.6g\n', WL_need);
        fprintf('  W_need      = %.6g m (%.3f um)\n', W_need, W_need*1e6);

        ts = datestr(now,'yyyymmdd_HHMMSS');
        outSizing = table( ...
            string(famName), string(key), string(L_sel), L_chosen, gm_target, normids_at_gm, Ids_need, WL_need, W_need, ...
            'VariableNames', {'family_autodetect','key','L_label','L_m','gmoverid_target','abs_normids_at_gm','Ids_need_A','W_over_L_need','W_need_m'} );
        outSizingCsv = fullfile(outputDir, sprintf('sizing_%s_%s.csv', lower(strrep(key,' ','_')), ts));
        writetable(outSizing, outSizingCsv);
        fprintf('Saved sizing row to: %s\n', outSizingCsv);

    catch ME
        fprintf(2, '\n[ERROR] %s\n\n', ME.message);
    end

    answ = lower(strtrim(input('\nDo another sizing run? (y/n): ', 's')));
    keepGoing = startsWith(answ, 'y');
    fprintf('\n');
end

fprintf('Exiting sizing loop.\n');

%% =========================================================================
%  Local functions
% =========================================================================
function res = processOneSet(inputDir, outputDir, name, prefix, L_labels_csv, L_values_csv, L_labels_legend)
    fn_gmoverid = fullfile(inputDir, sprintf('gmoverid_vs_vg_multL_%s.csv', name));
    fn_gmro     = fullfile(inputDir, sprintf('gmro_vs_vg_multL_%s.csv', name));
    fn_normids  = fullfile(inputDir, sprintf('normids_vs_vg_multL_%s.csv', name));

    T_gmoverid = readCsvAsNumericTable(fn_gmoverid);
    T_gmro     = readCsvAsNumericTable(fn_gmro);
    T_normids  = readCsvAsNumericTable(fn_normids);

    vg_ref = T_gmoverid{:,1};
    gmoverid_mat = T_gmoverid{:,2:end};

    gmro_mat    = interpToVgGrid(T_gmro,    vg_ref);
    normids_mat = interpToVgGrid(T_normids, vg_ref);

    % ABS normalized Ids everywhere
    normids_mat = abs(normids_mat);

    % Sanity check columns
    nL_found = size(gmoverid_mat, 2);
    nL_expected = numel(L_labels_csv);

    if nL_found ~= nL_expected
        warning('%s %s: expected %d L cols but found %d. Using min of both.', prefix, name, nL_expected, nL_found);
        nUse = min(nL_expected, nL_found);

        gmoverid_mat = gmoverid_mat(:, 1:nUse);
        gmro_mat     = gmro_mat(:,     1:nUse);
        normids_mat  = normids_mat(:,  1:nUse);

        L_labels_csv_use = L_labels_csv(1:nUse);
        L_values_csv_use = L_values_csv(1:nUse);
    else
        L_labels_csv_use = L_labels_csv;
        L_values_csv_use = L_values_csv;
        nUse = nL_found;
    end

    % Map desired legend order -> actual indices
    plotOrderIdx = [];
    plotOrderLbl = {};
    for k = 1:numel(L_labels_legend)
        lbl = L_labels_legend{k};
        idx = find(strcmp(L_labels_csv_use, lbl), 1, 'first');
        if ~isempty(idx)
            plotOrderIdx(end+1) = idx; %#ok<AGROW>
            plotOrderLbl{end+1} = lbl; %#ok<AGROW>
        end
    end
    if isempty(plotOrderIdx)
        warning('%s %s: legend mapping failed; falling back to CSV order.', prefix, name);
        plotOrderIdx = 1:nUse;
        plotOrderLbl = L_labels_csv_use;
    end

    plotOrderLm = zeros(numel(plotOrderIdx), 1);
    for ii = 1:numel(plotOrderIdx)
        plotOrderLm(ii) = L_values_csv_use(plotOrderIdx(ii));
    end

    % ----- CSV 1: gmoverid vs normids -----
    rows = table();
    for ii = 1:numel(plotOrderIdx)
        i = plotOrderIdx(ii);
        x = normids_mat(:,i); % abs
        y = gmoverid_mat(:,i);
        vg = vg_ref;
        Lm = L_values_csv_use(i);

        Ti = table( ...
            repmat(string(plotOrderLbl{ii}), numel(vg), 1), ...
            repmat(Lm,                      numel(vg), 1), ...
            vg, x, y, ...
            'VariableNames', {'L_label','L_m','Vg','normids','gmoverid'} );
        rows = [rows; Ti]; %#ok<AGROW>
    end

    outCsv1 = fullfile(outputDir, sprintf('%s_gmoverid_vs_normids_multL_%s.csv', prefix, name));
    writetable(rows, outCsv1);

    % ----- Plot 1: gmoverid vs normids (log x) -----
    fig1 = figure('Name', sprintf('%s_gmoverid_vs_normids_%s', prefix, name), 'Color', 'w');
    ax1 = axes(fig1);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

    h1 = gobjects(numel(plotOrderIdx), 1);
    for ii = 1:numel(plotOrderIdx)
        i = plotOrderIdx(ii);

        x = normids_mat(:,i);  x = x(:);
        y = gmoverid_mat(:,i); y = y(:);

        m = isfinite(x) & isfinite(y) & (x > 0);
        x = x(m); y = y(m);

        if numel(x) < 2
            h1(ii) = gobjects(1);
            continue;
        end

        [x, idx] = sort(x); y = y(idx);

        hTmp = plot(ax1, x, y, 'LineWidth', 1.2, 'DisplayName', plotOrderLbl{ii});
        if isempty(hTmp), h1(ii) = gobjects(1); else, h1(ii) = hTmp(1); end
    end

    set(ax1, 'XScale', 'log');
    xlabel(ax1, '|\it{I_D}\rm| normalized  (|I_D| / (W/L))');
    ylabel(ax1, 'g_m / I_D');
    title(ax1, sprintf('%s: g_m/I_D vs |Normalized I_D|  (%s)', prefix, upper(name)));

    good = isgraphics(h1);
    legend(ax1, h1(good), plotOrderLbl(good), 'Location', 'best');
    saveas(fig1, fullfile(outputDir, sprintf('%s_gmoverid_vs_normids_multL_%s.png', prefix, name)));

    % ----- CSV 2: gmro vs gmoverid -----
    rows2 = table();
    for ii = 1:numel(plotOrderIdx)
        i = plotOrderIdx(ii);
        x = gmoverid_mat(:,i);
        y = gmro_mat(:,i);
        vg = vg_ref;
        Lm = L_values_csv_use(i);

        Ti = table( ...
            repmat(string(plotOrderLbl{ii}), numel(vg), 1), ...
            repmat(Lm,                      numel(vg), 1), ...
            vg, x, y, ...
            'VariableNames', {'L_label','L_m','Vg','gmoverid','gmro'} );
        rows2 = [rows2; Ti]; %#ok<AGROW>
    end

    outCsv2 = fullfile(outputDir, sprintf('%s_gmro_vs_gmoverid_multL_%s.csv', prefix, name));
    writetable(rows2, outCsv2);

    % ----- Plot 2: gmro vs gmoverid -----
    fig2 = figure('Name', sprintf('%s_gmro_vs_gmoverid_%s', prefix, name), 'Color', 'w');
    ax2 = axes(fig2);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

    h2 = gobjects(numel(plotOrderIdx), 1);
    for ii = 1:numel(plotOrderIdx)
        i = plotOrderIdx(ii);

        x = gmoverid_mat(:,i); x = x(:);
        y = gmro_mat(:,i);     y = y(:);

        m = isfinite(x) & isfinite(y);
        x = x(m); y = y(m);

        if numel(x) < 2
            h2(ii) = gobjects(1);
            continue;
        end

        [x, idx] = sort(x); y = y(idx);

        hTmp = plot(ax2, x, y, 'LineWidth', 1.2, 'DisplayName', plotOrderLbl{ii});
        if isempty(hTmp), h2(ii) = gobjects(1); else, h2(ii) = hTmp(1); end
    end

    xlabel(ax2, 'g_m / I_D');
    ylabel(ax2, 'g_m r_o');
    title(ax2, sprintf('%s: g_m r_o vs g_m/I_D  (%s)', prefix, upper(name)));

    good2 = isgraphics(h2);
    legend(ax2, h2(good2), plotOrderLbl(good2), 'Location', 'best');
    saveas(fig2, fullfile(outputDir, sprintf('%s_gmro_vs_gmoverid_multL_%s.png', prefix, name)));

    fprintf('Done %s: %s\n  Wrote: %s\n  Wrote: %s\n', prefix, name, outCsv1, outCsv2);

    % Store for sizing (legend order)
    res = struct();
    res.L_labels = plotOrderLbl(:);
    res.L_m      = plotOrderLm(:);
    res.gmoverid = cell(numel(plotOrderIdx), 1);
    res.normids  = cell(numel(plotOrderIdx), 1);

    for ii = 1:numel(plotOrderIdx)
        i = plotOrderIdx(ii);
        res.gmoverid{ii} = gmoverid_mat(:,i);
        res.normids{ii}  = normids_mat(:,i); % abs
    end
end

function T = readCsvAsNumericTable(filename)
    if ~isfile(filename)
        error('Missing file: %s', filename);
    end

    try
        opts = detectImportOptions(filename, 'NumHeaderLines', 0);
        for k = 1:numel(opts.VariableTypes)
            opts.VariableTypes{k} = 'double';
        end
        T = readtable(filename, opts);

        isNum = varfun(@isnumeric, T, 'OutputFormat', 'uniform');
        if ~all(isNum)
            error('Non-numeric detected, fallback.');
        end
    catch
        M = readmatrix(filename);
        if isempty(M) || size(M,2) < 2
            error('File %s did not read as a numeric matrix with >=2 columns.', filename);
        end
        T = array2table(M);
    end
end

function data_mat = interpToVgGrid(T_in, vg_target)
    vg_in = T_in{:,1};
    data_in = T_in{:,2:end};

    if isequal(vg_in, vg_target)
        data_mat = data_in;
        return;
    end

    nCol = size(data_in,2);
    data_mat = nan(numel(vg_target), nCol);

    for c = 1:nCol
        y = data_in(:,c);
        m = isfinite(vg_in) & isfinite(y);

        if nnz(m) >= 2
            [vg_sorted, idx] = sort(vg_in(m));
            y_sorted = y(m);
            y_sorted = y_sorted(idx);
            data_mat(:,c) = interp1(vg_sorted, y_sorted, vg_target, 'linear', NaN);
        end
    end
end

function normids_at_gm = interpNormIdsAtGm(gm_vec, normids_vec, gm_target)
    gm = gm_vec(:);
    ni = normids_vec(:);

    m = isfinite(gm) & isfinite(ni);
    gm = gm(m);
    ni = ni(m);

    if numel(gm) < 2
        normids_at_gm = NaN;
        return;
    end

    [gm, idx] = sort(gm);
    ni = ni(idx);

    [gm_u, ia] = unique(gm, 'stable');
    ni_u = ni(ia);

    if numel(gm_u) < 2
        normids_at_gm = NaN;
        return;
    end

    normids_at_gm = interp1(gm_u, ni_u, gm_target, 'linear', NaN);
end
