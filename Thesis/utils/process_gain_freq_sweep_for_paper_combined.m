close all

%%

% READ output of keysight_rigol_gain_freqsweep.py over N iterations for averaging
%run_iter = 3;
run_iter=1;
%% ---- ITAIL=1nA
figure
subplot(3,1,1)
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_1nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

vamp=0.1;
vofs=-0.06;

freq=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    freq(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

freq=freq';
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

%-----
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

%--SIM
csvFileName = 'sim_csv\bodeplot_itail_1nA_V3_27C.csv'  % V1 1st stage, V2, 2nd, V3, VPOS_AFE
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim = csv_val(:,1);
gain_sim = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part

csvFileName = 'sim_csv\bodeplot_itail_1nA_V3_CL_27C.csv'
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim_cl = csv_val(:,1);
gain_sim_cl = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part
%-----

%semilogx(mean_freq, 20*log10(mean_gain), 'b', 'LineWidth', 2);
semilogx(freq_sim_cl,20*log10(gain_sim_cl), 'MarkerFaceColor', 'r');
hold on;
semilogx(freq_sim,20*log10(gain_sim), 'MarkerFaceColor', 'b');
% extrapolate cl gain
Cin=172e-15; Cf=40e-15; Cpar=5e-15;
% Add square data points using the scatter function
mean_gain_cl = (Cin/(Cf+Cin)).*mean_gain./(1+(mean_gain.*(Cf/(Cf+Cin))));
mean_gain_cl_itail_1nA = mean_gain_cl; mean_gain_itail_1nA = mean_gain; % save for later combined plots
scatter(mean_freq, 20*log10(mean_gain_cl), 100, 's', 'filled', 'MarkerFaceColor', 'b'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
% Add square data points using the scatter function
scatter(mean_freq, 20*log10(mean_gain), 100, 's', 'filled', 'MarkerFaceColor', 'r'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
xlim([10e-6,10e3])
ylim([0.1*min(20*log10(mean_gain)),1.25*max(20*log10(mean_gain))])
%legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
%xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title(['a) I_{TAIL}=1nA']);

hold off;
micasplot
%% ---- ITAIL=10nA
subplot(3,1,2)
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_10nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

vamp=0.1;
vofs=-0.02;

freq=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    freq(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

freq=freq';
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

%-----
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

%--SIM
csvFileName = 'sim_csv\bodeplot_itail_10nA_V3_27C.csv'  % V1 1st stage, V2, 2nd, V3, VPOS_AFE
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim = csv_val(:,1);
gain_sim = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part

csvFileName = 'sim_csv\bodeplot_itail_10nA_V3_CL_27C.csv'
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim_cl = csv_val(:,1);
gain_sim_cl = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part
%-----
grid on;
%semilogx(mean_freq, 20*log10(mean_gain), 'b', 'LineWidth', 2);
semilogx(freq_sim_cl,20*log10(gain_sim_cl), 'MarkerFaceColor', 'r');
hold on;
semilogx(freq_sim,20*log10(gain_sim), 'MarkerFaceColor', 'b');
% extrapolate cl gain
Cin=172e-15; Cf=40e-15; Cpar=5e-15;
% Add square data points using the scatter function
mean_gain_cl = (Cin/(Cf+Cin)).*mean_gain./(1+(mean_gain.*(Cf/(Cf+Cin))));
mean_gain_cl_itail_10nA = mean_gain_cl; mean_gain_itail_10nA = mean_gain; % save for later combined plots
scatter(mean_freq, 20*log10(mean_gain_cl), 100, 's', 'filled', 'MarkerFaceColor', 'b'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
% Add square data points using the scatter function
scatter(mean_freq, 20*log10(mean_gain), 100, 's', 'filled', 'MarkerFaceColor', 'r'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
xlim([10e-6,10e3])
ylim([0.1*min(20*log10(mean_gain)),1.25*max(20*log10(mean_gain))])
%legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
%xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title(['b) I_{TAIL}=10nA']);


hold off;
micasplot


%% ---- ITAIL=50nA
subplot(3,1,3)
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_50nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

vamp=0.1;
vofs=-0.15;

freq=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    freq(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

freq=freq';
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

%-----
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

%--SIM
csvFileName = 'sim_csv\bodeplot_itail_50nA_V3_27C.csv'  % V1 1st stage, V2, 2nd, V3, VPOS_AFE
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim = csv_val(:,1);
gain_sim = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part

csvFileName = 'sim_csv\bodeplot_itail_50nA_V3_CL_27C.csv'
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim_cl = csv_val(:,1);
gain_sim_cl = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part
%-----

%semilogx(mean_freq, 20*log10(mean_gain), 'b', 'LineWidth', 2);
semilogx(freq_sim_cl,20*log10(gain_sim_cl), 'MarkerFaceColor', 'r');
hold on;
semilogx(freq_sim,20*log10(gain_sim), 'MarkerFaceColor', 'b');
% extrapolate cl gain
Cin=172e-15; Cf=40e-15; Cpar=5e-15;
% Add square data points using the scatter function
mean_gain_cl = (Cin/(Cf+Cin)).*mean_gain./(1+(mean_gain.*(Cf/(Cf+Cin))));
mean_gain_cl_itail_50nA = mean_gain_cl; mean_gain_itail_50nA = mean_gain; % save for later combined plots
scatter(mean_freq, 20*log10(mean_gain_cl), 100, 's', 'filled', 'MarkerFaceColor', 'b'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
% Add square data points using the scatter function
scatter(mean_freq, 20*log10(mean_gain), 100, 's', 'filled', 'MarkerFaceColor', 'r'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
xlim([10e-6,10e3])
ylim([0.1*min(20*log10(mean_gain)),1.25*max(20*log10(mean_gain))])
%legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title(['c) I_{TAIL}=50nA']);

hold off;
micasplot

hlegend = legend(gca,'Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL}','Measured A_{OL}','location','Best')
set(hlegend, 'FontSize', 10)


% Set the width and height of the figure (in pixels)
width = 1000;    % Specify your desired width here
height = 800;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,['output\00_forpaper_bodeplot_acrossitail.fig'])
saveas(gcf,['output\00_forpaper_bodeplot_acrossitail.png'])
%%
clear all

