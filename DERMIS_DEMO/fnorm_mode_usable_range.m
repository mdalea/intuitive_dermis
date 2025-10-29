% MATLAB script: separate figures for Cf=10fF and Cf=100fF
close all
clear; clc;

VTX = 3;                        % Transmit voltage [V]
Cf_list = [10e-15, 100e-15];    % Cf values
Cf_labels = ["Cf=10fF", "Cf=100fF"];
Cf_filenames = ["TBCAS25_fnorm_mode_usable_range_Cf10.png", ...
                "TBCAS25_fnorm_mode_usable_range_Cf100.png"];

% Sweep values for Ctot = Csp + Csn (assume Csp = Csn, so Ctot = 2*Csp)
Ctot = logspace(-18, -12, 800);   % ~0.01 aF to 1 pF
Csp  = Ctot/2;
Csn  = Ctot/2;

% Integrated sensor max (vertical marker)
Ctot_max_int = 6.138e-15 / 2;     % corrected value

for k = 1:2
    Cf = Cf_list(k);

    % Denominator
    den = Cf * (Ctot + 2*Cf);

    % Contributions to VOUT
    VoutA = VTX .* (2 .* Csp .* Csn) ./ den;       
    VoutB = VTX .* (Cf .* (Ctot)) ./ den;     
    Vout  = VoutA + VoutB;                         

    % ---- New figure ----
    figure('Color','w');

    % Main plot (log-log)
    loglog(Ctot, VoutA, 'LineWidth', 2); hold on;
    loglog(Ctot, VoutB, 'LineWidth', 2);
    loglog(Ctot, Vout,  'LineWidth', 2, 'LineStyle', ':');

    % Axes labels & limits
    grid on;
    xlabel('Csp + Csn [F]');
    ylabel(sprintf('VOUT (V) @ %s', Cf_labels(k)));
    ylim([0 20]);                                  % set limits first

    % Vertical thick dashed line at integrated sensor max
    xline(Ctot_max_int,  '--r', 'Max Csp+Csn', 'LineWidth',2.5, 'FontWeight', 'bold');

    % Horizontal line at ADC full-scale (differential)
    yline(6, '--r', 'ADC FS Vpp (differential)', ...
          'LineWidth', 2.5, 'LabelHorizontalAlignment', 'left', ...
          'LabelVerticalAlignment', 'bottom', 'FontWeight', 'bold');

    % Place "Max Csp+Csn" text at vertical middle, very close to the line
    yl = ylim;
    y_mid = mean(yl);
    x_text = Ctot_max_int * 1.03;                  % small log-scale offset to the right
    % text(x_text, y_mid, 'Max Csp+Csn', ...
    %      'FontWeight','bold', 'HorizontalAlignment','left', ...
    %      'VerticalAlignment','middle');

    legend( ...
        'VTX·2·Csp·Csn / (Cf·(Csp+Csn+2Cf))', ...
        'VTX·Cf·(Csp+Csn) / (Cf·(Csp+Csn+2Cf))', ...
        'Total VOUT', ...
        'Location','southwest');

    % ---- Inset axes: Relative error if Term A neglected ----
    err_rel_pct = (VoutA ./ Vout) * 100;
    err_at_max = interp1(Ctot, err_rel_pct, Ctot_max_int, 'linear');

    % Inset in lower-right corner
    ax = gca;
    pos = ax.Position;              
    inset_width  = 0.18 * pos(3);   
    inset_height = 0.55 * pos(4) * 0.9;   
    margin_x = 0.05 * pos(3);
    margin_y = 0.05 * pos(4);
    inset_x = pos(1) + pos(3) - inset_width - margin_x;
    inset_y = pos(2) + margin_y + 0.20*inset_height;

    ax_in = axes('Position', [inset_x, inset_y, inset_width, inset_height]);

    semilogx(Ctot, err_rel_pct, 'LineWidth', 1.5); hold on;
    xline(Ctot_max_int, '--k', 'LineWidth', 1.5);
    grid on; box on;
    ylim([0 20]);

    % Labels inside inset
    text(0.02, 0.95, 'Error (%)', 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 8);
    % Error value at middle-left
    text(0.05, 0.50, sprintf('%.2f%%', err_at_max), 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
         'FontSize', 8, 'FontAngle','italic');

    set(ax_in, 'FontSize', 8);
    xlim([min(Ctot) max(Ctot)]);

    % Optional style tweak (if you use micasplot)
    micasplot
    ax_all = findall(gcf,'Type','axes');  for axx = ax_all', axx.FontSize = axx.FontSize * 0.95; end
    leg_all = findall(gcf,'Type','legend'); for lg = leg_all', lg.FontSize = lg.FontSize * 0.95; end
    txt_all = findall(gcf,'Type','text');  for tx = txt_all', tx.FontSize = tx.FontSize * 0.95; end

    % Save figure
    width = 800; height = 300;
    set(gcf, 'Position', [100, 100, width, height]);
    saveas(gcf, Cf_filenames(k));
end

