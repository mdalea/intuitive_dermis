% Define capacitances and transmitter voltage
Csp = 1.7e-12;  %6.12e-15; %1.7e-12;  % Capacitance Csp in Farads
Csn = 1.7e-12;  %6.12e-15; %1.7e-12;  % Capacitance Csn in Farads
Cf = 40e-15;    % Capacitance Cf in Farads
Catt = 1.7e-12; % 0; %1.7e-12; % Capacitance Catt in Farads
Vtx = 3;        % Transmitter voltage

% Configure measurement parameters
samplingRate = 500e3; % Sampling rate (500 kHz in this example)
windowTime = 2;
windowSize = samplingRate * windowTime;
overlap = windowSize / 2;
N = 3; % Number of files to average

% Initialize accumulators for PSD estimates
psdEstimate_dt_accum = 0;
psdEstimate_ct_accum = 0;

for i = 0:(N-1) % Loop through the N CSV files
    baseName = 'sample_3_vtune_1p15_vcmin_1p513_vcm_1p5_vddafe_3p0_ibias_40n_fixed_nosig_noisepsd_dt'
    baseName = 'sample_4_ibias_50n_vtune_0p2_fchop_2k_nosig2_noisepsd_dt'

    % ----- Get AFE output (DT mode)
    csvFileName_voutp_dt = sprintf('output\\ID_%s__wav_voutp_x0-%d.csv', baseName,i)
    csvFileName_voutn_dt = sprintf('output\\ID_%s__wav_voutn_x0-%d.csv', baseName,i)

    waveform_voutp_dt = csvread(csvFileName_voutp_dt);
    waveform_voutn_dt = csvread(csvFileName_voutn_dt);
    waveform_dt = waveform_voutp_dt - waveform_voutn_dt;
    
    % Calculate PSD for DT mode
    [psdEstimate_dt, freq_dt] = pwelch(waveform_dt(:,2), windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided');

    baseName = 'sample_3_vtune_0p35_vcmin_1p513_vcm_1p5_vddafe_3p0_ibias_40n_fixed_nosig_noisepsd_ct'
    baseName = 'sample_4_ibias_50n_vtune_0p0_fchop_2k_nosig3_noisepsd_ct'
    
    % ----- Get AFE output (CT mode)
    csvFileName_voutp_ct = sprintf('output\\ID_%s__wav_voutp_x0-%d.csv', baseName,i)
    csvFileName_voutn_ct = sprintf('output\\ID_%s__wav_voutn_x0-%d.csv', baseName,i)

    waveform_voutp_ct = csvread(csvFileName_voutp_ct);
    waveform_voutn_ct = csvread(csvFileName_voutn_ct);
    waveform_ct = waveform_voutp_ct - waveform_voutn_ct;
    
    % Calculate PSD for CT mode
    [psdEstimate_ct, freq_ct] = pwelch(waveform_ct(:,2), windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided');
    
    % Accumulate PSD estimates
    psdEstimate_dt_accum = psdEstimate_dt_accum + psdEstimate_dt;
    psdEstimate_ct_accum = psdEstimate_ct_accum + psdEstimate_ct;
end

% Average the PSD estimates
psdEstimate_dt_avg = psdEstimate_dt_accum / N;
psdEstimate_ct_avg = psdEstimate_ct_accum / N;

% Calculate total capacitance
totalCap = (Csp + Csn + 2*Cf + 2*Catt);

%% ------- Capacitance-Referred Input PSD Calculation

% Capacitance-referred input PSD (Cs = [Vout_PSD * (Csp + Csn + 2*Cf + 2*Catt)] / Vtx)
psdEstimate_dt_Cin = (psdEstimate_dt_avg * totalCap^2) / Vtx^2;  % DT mode
psdEstimate_ct_Cin = (psdEstimate_ct_avg * totalCap^2) / Vtx^2;  % CT mode

% Convert capacitance-referred PSD to F/sqrt(Hz)
psdEstimate_dt_Cin_F = sqrt(psdEstimate_dt_Cin);  % DT mode
psdEstimate_ct_Cin_F = sqrt(psdEstimate_ct_Cin);  % CT mode

% Find indices corresponding to the frequency range of interest
flo = 1; % Specify the lower frequency of the bandwidth (in Hz)
fhi = 300; %1e3; % Specify the higher frequency of the bandwidth (in Hz)
idxRange_dt = find(freq_dt >= flo & freq_dt <= fhi);
idxRange_ct = find(freq_ct >= flo & freq_ct <= fhi);

% Integrate the PSD over the specified frequency range
powerInBand_dt = sum(psdEstimate_dt_Cin(idxRange_dt));
powerInBand_ct = sum(psdEstimate_ct_Cin(idxRange_ct));

% Calculate the RMS capacitance based on the power
rmsCap_dt = sqrt(powerInBand_dt) ;
rmsCap_ct = sqrt(powerInBand_ct) ;

%% ------- Plotting Capacitance-Referred Input PSD (F^2/Hz)

figure;
loglog(freq_dt, psdEstimate_dt_Cin, 'LineWidth', 1.5); xlim([1 30e3])
hold on;
loglog(freq_ct, psdEstimate_ct_Cin, 'LineWidth', 1.5); xlim([1 30e3])
xlabel('Frequency (Hz)');
ylabel('Capacitance-Referred PSD (F^2/Hz)');
legend('Chopper on', 'Chopper off/ IDLE mode', 'Location', 'southwest')
title('Capacitance-Referred Input PSD (F^2/Hz)');
grid on;
hold off;

% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0005*freq_dt(round(length(freq_dt) / 2));
y_pos = min(psdEstimate_ct_Cin * 2e7);
text(x_pos, y_pos, sprintf('(BW: %d Hz to %d Hz) \nChopper on: F_{noise RMS} = %.3f fF \nChopper off: F_{noise RMS} = %.3f fF', flo, fhi, 1e15*rmsCap_dt, 1e15*rmsCap_ct), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');


micasplot
% Save the figure
width = 800;
height = 600;
set(gcf, 'Position', [100, 100, width, height]);
% saveas(gcf, [csvFileName_voutp_dt,'capacitance_referred_psd.fig']);
% saveas(gcf, [csvFileName_voutp_dt,'capacitance_referred_psd.png']);
saveas(gcf, ['output\',baseName,'_capacitance_referred_psd.fig']);
saveas(gcf, ['output\',baseName,'_capacitance_referred_psd.png']);

%% ------- Plotting Capacitance-Referred Input PSD (F/sqrt(Hz))

figure;
loglog(freq_dt, psdEstimate_dt_Cin_F, 'LineWidth', 1.5); xlim([1 30e3])
hold on;
loglog(freq_ct, psdEstimate_ct_Cin_F, 'LineWidth', 1.5); xlim([1 30e3])
xlabel('Frequency (Hz)');
ylabel('Capacitance-Referred PSD (F/sqrt(Hz))');
legend('Chopper on', 'Chopper off/ IDLE mode', 'Location', 'southwest')
title('Capacitance-Referred Input PSD (F/sqrt(Hz))');
grid on;
hold off;


% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0005*freq_dt(round(length(freq_dt) / 2));
y_pos = min(psdEstimate_ct_Cin_F * 2e3);
text(x_pos, y_pos, sprintf('(BW: %d Hz to %d Hz) \nChopper on: F_{noise RMS} = %.3f fF \nChopper off: F_{noise RMS} = %.3f fF', flo, fhi, 1e15*rmsCap_dt, 1e15*rmsCap_ct), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

micasplot
% Save the figure
set(gcf, 'Position', [100, 100, width, height]);
% saveas(gcf, [csvFileName_voutp_dt,'capacitance_referred_psd_F.fig']);
% saveas(gcf, [csvFileName_voutp_dt,'capacitance_referred_psd_F.png']);
saveas(gcf, ['output\',baseName,'_capacitance_referred_psd_F.fig']);
saveas(gcf, ['output\',baseName,'_capacitance_referred_psd_F.png']);
