close all

%%

% READ output of keysight_rigol_gain_freqsweep.py over N iterations for averaging
%run_iter = 3;
run_iter=1;
%% ---- ITAIL=1nA
%run_id = 'sample_I2_covered_for_paper_itail_1nA'
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_1nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

%used in TSENSE
fin_arr = logspace(0,4,21); % 5 points/dec over 4 decade
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



figure
loglog(freq,gain)
ylim([0.5*min(min(gain)),1.5*max(max(gain))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])

%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get an estimate of the passband gain
gain_max = max(mean_gain);

% Plot the means as a solid line
loglog(mean_freq, mean_gain, 'b', 'LineWidth', 2);
hold on;
ylim([0.5*min(mean_gain),1.5*max(mean_gain)])

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','southwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 10*max(mean_gain);
%text(x_pos, y_pos, sprintf('Measured Passband Gain = %.2f V/V', gain_max), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

micasplot
hold off;
% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])

fig_itail_1nA = gcf;
ax_itail_1nA = gca;


%------ plot without error bars (for paper)
figure

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
ylim([0.1*min(20*log10(mean_gain)),2*max(20*log10(mean_gain))])
legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);


% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 1.25*max(20*log10(mean_gain));
%text(x_pos, y_pos, sprintf('Passband Gain = %.2f dB', 20*log10(gain_max)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-points','.fig'])
saveas(gcf,[baseName,'multitrials-points','.png'])
% plot vout and vin -> try and see if vin is bandlimited by eval-ADA4098
figure
loglog(freq,vout_pp)
hold on;
loglog(freq,vin_pp)
ylim([0.5*min(min(vin_pp)),1.5*max(max(vout_pp))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Output p-p Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin','.png'])
hold off;
%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_vout_pp = mean(vout_pp,2); % get average across iterations 2nd dim
std_vout_pp = std(vout_pp,0,2); % get sigma
mean_vin_pp = mean(vin_pp,2); % get average across iterations 2nd dim
std_vin_pp = std(vin_pp,0,2); % get sigma

% Plot the means as a solid line
loglog(mean_freq, mean_vout_pp, 'b', 'LineWidth', 2);
hold on;
loglog(mean_freq, mean_vin_pp, 'b', 'LineWidth', 2);
ylim([0.5*min(mean_vin_pp),1.5*max(mean_vout_pp)])



% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_vout_pp, std_vout_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(mean_freq, mean_vin_pp, std_vin_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean (Vout,pp)','Mean (Vin,pp)', 'Standard Deviation (Vout,pp)', 'Standard Deviation (Vin,pp)','location','southwest');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.png'])


% Plot temperature over iterations and sweep values
baseName_1 = ['output\ID_', run_id, '_bodeplot__tsense_ave_fin_']
baseName_2 = ['_vpp_',num2str(vamp),'_vofs_', num2str(vofs)]


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(fin_arr)
        if (mod(fin_arr(j),1) == 0) % integer
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.1f'),baseName_2,'-',num2str(i-1),'.csv'] % that pesky case of num2str conversion in python
        else
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.6f'),baseName_2,'-',num2str(i-1),'.csv']            
        end
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
fin_arr=fin_arr'; 
figure
hold on;
plot(fin_arr,tsense_arr)
set(gca,'Xscale','log') % yeah, that is indeed easy!
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('TSENSE output (V)');
title(['TSENSE across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5');

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense','.fig'])
saveas(gcf,[baseName,'multitrials-tsense','.png'])
% See general trend of TSENSE

% Calculate the mean and standard deviation along dimension 2 (columns) of tsense_arr
means = mean(tsense_arr, 1);
sigmas = std(tsense_arr, 0, 1);

% Create the figure
figure;

% Plot the means as a solid line against dimension 2
plot(1:size(tsense_arr, 2), means, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation) against dimension 2
errorbar(1:size(tsense_arr, 2), means, sigmas, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Iterations');
ylabel('TSENSE (V)');
title('Mean and Standard Deviation of TSENSE across run iterations');

% Add legend
legend('Mean', 'Standard Deviation');

% Limit the X-axis ticks to show only the values of dimension 2
xticks(1:size(tsense_arr, 2));

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense-trend','.fig'])
saveas(gcf,[baseName,'multitrials-tsense-trend','.png'])

%% ---- ITAIL=10nA
run_id = 'sample_I2_covered_for_paper_itail_10nA'
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_10nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

%used in TSENSE
fin_arr = logspace(0,4,21); % 5 points/dec over 4 decade
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



figure
loglog(freq,gain)
ylim([0.5*min(min(gain)),1.5*max(max(gain))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])

%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get an estimate of the passband gain
gain_max = max(mean_gain);

% Plot the means as a solid line
loglog(mean_freq, mean_gain, 'b', 'LineWidth', 2);
hold on;
ylim([0.5*min(mean_gain),1.5*max(mean_gain)])

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','southwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 10*max(mean_gain);
%text(x_pos, y_pos, sprintf('Measured Passband Gain = %.2f V/V', gain_max), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

micasplot
hold off;
% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])

fig_itail_10nA = gcf;
ax_itail_10nA = gca;
%------ plot without error bars (for paper)
figure

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
ylim([0.1*min(20*log10(mean_gain)),2*max(20*log10(mean_gain))])
legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);


% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 1.25*max(20*log10(mean_gain));
%text(x_pos, y_pos, sprintf('Passband Gain = %.2f dB', 20*log10(gain_max)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-points','.fig'])
saveas(gcf,[baseName,'multitrials-points','.png'])
% plot vout and vin -> try and see if vin is bandlimited by eval-ADA4098
figure
loglog(freq,vout_pp)
hold on;
loglog(freq,vin_pp)
ylim([0.5*min(min(vin_pp)),1.5*max(max(vout_pp))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Output p-p Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin','.png'])
hold off;
%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_vout_pp = mean(vout_pp,2); % get average across iterations 2nd dim
std_vout_pp = std(vout_pp,0,2); % get sigma
mean_vin_pp = mean(vin_pp,2); % get average across iterations 2nd dim
std_vin_pp = std(vin_pp,0,2); % get sigma

% Plot the means as a solid line
loglog(mean_freq, mean_vout_pp, 'b', 'LineWidth', 2);
hold on;
loglog(mean_freq, mean_vin_pp, 'b', 'LineWidth', 2);
ylim([0.5*min(mean_vin_pp),1.5*max(mean_vout_pp)])



% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_vout_pp, std_vout_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(mean_freq, mean_vin_pp, std_vin_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean (Vout,pp)','Mean (Vin,pp)', 'Standard Deviation (Vout,pp)', 'Standard Deviation (Vin,pp)','location','southwest');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.png'])

% Plot temperature over iterations and sweep values
baseName_1 = ['output\ID_', run_id, '_bodeplot__tsense_ave_fin_']
baseName_2 = ['_vpp_',num2str(vamp),'_vofs_', num2str(vofs)]


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(fin_arr)
        if (mod(fin_arr(j),1) == 0) % integer
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.1f'),baseName_2,'-',num2str(i-1),'.csv'] % that pesky case of num2str conversion in python
        else
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.6f'),baseName_2,'-',num2str(i-1),'.csv']            
        end
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
fin_arr=fin_arr'; 
figure
hold on;
plot(fin_arr,tsense_arr)
set(gca,'Xscale','log') % yeah, that is indeed easy!
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('TSENSE output (V)');
title(['TSENSE across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5');

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense','.fig'])
saveas(gcf,[baseName,'multitrials-tsense','.png'])
% See general trend of TSENSE

% Calculate the mean and standard deviation along dimension 2 (columns) of tsense_arr
means = mean(tsense_arr, 1);
sigmas = std(tsense_arr, 0, 1);

% Create the figure
figure;

% Plot the means as a solid line against dimension 2
plot(1:size(tsense_arr, 2), means, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation) against dimension 2
errorbar(1:size(tsense_arr, 2), means, sigmas, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Iterations');
ylabel('TSENSE (V)');
title('Mean and Standard Deviation of TSENSE across run iterations');

% Add legend
legend('Mean', 'Standard Deviation');

% Limit the X-axis ticks to show only the values of dimension 2
xticks(1:size(tsense_arr, 2));

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense-trend','.fig'])
saveas(gcf,[baseName,'multitrials-tsense-trend','.png'])

%% ---- ITAIL=100nA
run_id = 'sample_I2_covered_for_paper_itail_100nA'
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_100nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

%used in TSENSE
fin_arr = logspace(0,4,21); % 5 points/dec over 4 decade
vamp=0.1;
vofs=-0.12;

freq=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ]
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



figure
loglog(freq,gain)
ylim([0.5*min(min(gain)),1.5*max(max(gain))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])

%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get an estimate of the passband gain
gain_max = max(mean_gain);

% Plot the means as a solid line
loglog(mean_freq, mean_gain, 'b', 'LineWidth', 2);
hold on;
ylim([0.5*min(mean_gain),1.5*max(mean_gain)])

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','southwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 10*max(mean_gain);
%text(x_pos, y_pos, sprintf('Measured Passband Gain = %.2f V/V', gain_max), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

micasplot
hold off;
% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])

fig_itail_100nA = gcf;
ax_itail_100nA = gca;
%------ plot without error bars (for paper)
figure

%--SIM
csvFileName = 'sim_csv\bodeplot_itail_100nA_V3_27C.csv'  % V1 1st stage, V2, 2nd, V3, VPOS_AFE
csv_val = csvread(csvFileName,1,0); % skip first line of names
freq_sim = csv_val(:,1);
gain_sim = sqrt(csv_val(:,2).^2 + csv_val(:,3).^2); % negative Real part

csvFileName = 'sim_csv\bodeplot_itail_100nA_V3_CL_27C.csv'
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
mean_gain_cl_itail_100nA = mean_gain_cl; mean_gain_itail_100nA = mean_gain; % save for later combined plots
scatter(mean_freq, 20*log10(mean_gain_cl), 100, 's', 'filled', 'MarkerFaceColor', 'b'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
% Add square data points using the scatter function
scatter(mean_freq, 20*log10(mean_gain), 100, 's', 'filled', 'MarkerFaceColor', 'r'); % 's' specifies square markers, 'filled' fills the squares, 'r' is the marker face color (red)
xlim([10e-6,10e3])
ylim([0.1*min(20*log10(mean_gain)),2*max(20*log10(mean_gain))])
legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);


% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 1.25*max(20*log10(mean_gain));
%text(x_pos, y_pos, sprintf('Passband Gain = %.2f dB', 20*log10(gain_max)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-points','.fig'])
saveas(gcf,[baseName,'multitrials-points','.png'])
% plot vout and vin -> try and see if vin is bandlimited by eval-ADA4098
figure
loglog(freq,vout_pp)
hold on;
loglog(freq,vin_pp)
ylim([0.5*min(min(vin_pp)),1.5*max(max(vout_pp))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Output p-p Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin','.png'])
hold off;
%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_vout_pp = mean(vout_pp,2); % get average across iterations 2nd dim
std_vout_pp = std(vout_pp,0,2); % get sigma
mean_vin_pp = mean(vin_pp,2); % get average across iterations 2nd dim
std_vin_pp = std(vin_pp,0,2); % get sigma

% Plot the means as a solid line
loglog(mean_freq, mean_vout_pp, 'b', 'LineWidth', 2);
hold on;
loglog(mean_freq, mean_vin_pp, 'b', 'LineWidth', 2);
ylim([0.5*min(mean_vin_pp),1.5*max(mean_vout_pp)])



% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_vout_pp, std_vout_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(mean_freq, mean_vin_pp, std_vin_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean (Vout,pp)','Mean (Vin,pp)', 'Standard Deviation (Vout,pp)', 'Standard Deviation (Vin,pp)','location','southwest');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.png'])

% Plot temperature over iterations and sweep values
baseName_1 = ['output\ID_', run_id, '_bodeplot__tsense_ave_fin_']
baseName_2 = ['_vpp_',num2str(vamp),'_vofs_', num2str(vofs)]


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(fin_arr)
        if (mod(fin_arr(j),1) == 0) % integer
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.1f'),baseName_2,'-',num2str(i-1),'.csv'] % that pesky case of num2str conversion in python
        else
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.6f'),baseName_2,'-',num2str(i-1),'.csv']            
        end
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
fin_arr=fin_arr'; 
figure
hold on;
plot(fin_arr,tsense_arr)
set(gca,'Xscale','log') % yeah, that is indeed easy!
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('TSENSE output (V)');
title(['TSENSE across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5');

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense','.fig'])
saveas(gcf,[baseName,'multitrials-tsense','.png'])
% See general trend of TSENSE

% Calculate the mean and standard deviation along dimension 2 (columns) of tsense_arr
means = mean(tsense_arr, 1);
sigmas = std(tsense_arr, 0, 1);

% Create the figure
figure;

% Plot the means as a solid line against dimension 2
plot(1:size(tsense_arr, 2), means, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation) against dimension 2
errorbar(1:size(tsense_arr, 2), means, sigmas, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Iterations');
ylabel('TSENSE (V)');
title('Mean and Standard Deviation of TSENSE across run iterations');

% Add legend
legend('Mean', 'Standard Deviation');

% Limit the X-axis ticks to show only the values of dimension 2
xticks(1:size(tsense_arr, 2));

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense-trend','.fig'])
saveas(gcf,[baseName,'multitrials-tsense-trend','.png'])

%% ---- ITAIL=50nA
run_id = 'sample_I2_covered_for_paper_itail_50nA'
run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_50nA'
baseName = ['output\ID_', run_id, '_bodeplot__gain_vs_fin_sweep-']

%used in TSENSE
fin_arr = logspace(0,4,21); % 5 points/dec over 4 decade
vamp=0.1;
vofs=-0.15;

freq=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ]
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



figure
loglog(freq,gain)
ylim([0.5*min(min(gain)),1.5*max(max(gain))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])

%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get an estimate of the passband gain
gain_max = max(mean_gain);

% Plot the means as a solid line
loglog(mean_freq, mean_gain, 'b', 'LineWidth', 2);
hold on;
ylim([0.5*min(mean_gain),1.5*max(mean_gain)])

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','southwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 10*max(mean_gain);
%text(x_pos, y_pos, sprintf('Measured Passband Gain = %.2f V/V', gain_max), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

micasplot
hold off;
% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])

fig_itail_50nA = gcf;
ax_itail_50nA = gca;
%------ plot without error bars (for paper)
figure

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
ylim([0.1*min(20*log10(mean_gain)),2*max(20*log10(mean_gain))])
legend('Simulated A_{CL}','Simulated A_{OL}','Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);


% Add text annotation for vofs_opt at the bottom center
x_pos = 0.0001*mean_freq(round(length(mean_freq) / 2));
y_pos = 1.25*max(20*log10(mean_gain));
%text(x_pos, y_pos, sprintf('Passband Gain = %.2f dB', 20*log10(gain_max)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);

saveas(gcf,[baseName,'multitrials-points','.fig'])
saveas(gcf,[baseName,'multitrials-points','.png'])
% plot vout and vin -> try and see if vin is bandlimited by eval-ADA4098
figure
loglog(freq,vout_pp)
hold on;
loglog(freq,vin_pp)
ylim([0.5*min(min(vin_pp)),1.5*max(max(vout_pp))])

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Output p-p Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin','.png'])
hold off;
%-----
figure
mean_freq = mean(freq,2);  % get mean across dim=2 just to have it 1 column; 
mean_vout_pp = mean(vout_pp,2); % get average across iterations 2nd dim
std_vout_pp = std(vout_pp,0,2); % get sigma
mean_vin_pp = mean(vin_pp,2); % get average across iterations 2nd dim
std_vin_pp = std(vin_pp,0,2); % get sigma

% Plot the means as a solid line
loglog(mean_freq, mean_vout_pp, 'b', 'LineWidth', 2);
hold on;
loglog(mean_freq, mean_vin_pp, 'b', 'LineWidth', 2);
ylim([0.5*min(mean_vin_pp),1.5*max(mean_vout_pp)])



% Plot the sigmas as error bars (standard deviation)
errorbar(mean_freq, mean_vout_pp, std_vout_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(mean_freq, mean_vin_pp, std_vin_pp, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Vpp (V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean (Vout,pp)','Mean (Vin,pp)', 'Standard Deviation (Vout,pp)', 'Standard Deviation (Vin,pp)','location','southwest');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-vout-vin-boxplot','.png'])

% Plot temperature over iterations and sweep values
baseName_1 = ['output\ID_', run_id, '_bodeplot__tsense_ave_fin_']
baseName_2 = ['_vpp_',num2str(vamp),'_vofs_', num2str(vofs)]


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(fin_arr)
        if (mod(fin_arr(j),1) == 0) % integer
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.1f'),baseName_2,'-',num2str(i-1),'.csv'] % that pesky case of num2str conversion in python
        else
            csvFileName = [baseName_1,num2str(fin_arr(j),'%.6f'),baseName_2,'-',num2str(i-1),'.csv']            
        end
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
fin_arr=fin_arr'; 
figure
hold on;
plot(fin_arr,tsense_arr)
set(gca,'Xscale','log') % yeah, that is indeed easy!
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('TSENSE output (V)');
title(['TSENSE across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5');

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense','.fig'])
saveas(gcf,[baseName,'multitrials-tsense','.png'])
% See general trend of TSENSE

% Calculate the mean and standard deviation along dimension 2 (columns) of tsense_arr
means = mean(tsense_arr, 1);
sigmas = std(tsense_arr, 0, 1);

% Create the figure
figure;

% Plot the means as a solid line against dimension 2
plot(1:size(tsense_arr, 2), means, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation) against dimension 2
errorbar(1:size(tsense_arr, 2), means, sigmas, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Iterations');
ylabel('TSENSE (V)');
title('Mean and Standard Deviation of TSENSE across run iterations');

% Add legend
legend('Mean', 'Standard Deviation');

% Limit the X-axis ticks to show only the values of dimension 2
xticks(1:size(tsense_arr, 2));

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-tsense-trend','.fig'])
saveas(gcf,[baseName,'multitrials-tsense-trend','.png'])


%%
%%COMBINED ACROSS ITAILS
figure
% Define the new colors using the RGB triplets
color_itail_1nA = [0 0.4470 0.7410];   % Dark blue
color_itail_10nA = [0.8500 0.3250 0.0980];  % Dark orange
color_itail_50nA = [0.9290 0.6940 0.1250]; % Dark yellow

% Scatter plots with itail_1nA data (set to blue)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_1nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_1nA); hold on;
scatter(mean_freq, 20*log10(mean_gain_itail_1nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_1nA);
% Tracing line for itail_1nA data
plot(mean_freq, 20*log10(mean_gain_cl_itail_1nA), 'Color', color_itail_1nA);
plot(mean_freq, 20*log10(mean_gain_itail_1nA), 'Color', color_itail_1nA);

% Scatter plots with itail_10nA data (set to green)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_10nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_10nA);
scatter(mean_freq, 20*log10(mean_gain_itail_10nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_10nA);
% Tracing line for itail_10nA data
plot(mean_freq, 20*log10(mean_gain_cl_itail_10nA), 'Color', color_itail_10nA);
plot(mean_freq, 20*log10(mean_gain_itail_10nA), 'Color', color_itail_10nA);

% Scatter plots with itail_50nA data (set to red)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_50nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_50nA);
scatter(mean_freq, 20*log10(mean_gain_itail_50nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_50nA);
% Tracing line for itail_50nA data
plot(mean_freq, 20*log10(mean_gain_cl_itail_50nA), 'Color', color_itail_50nA);
plot(mean_freq, 20*log10(mean_gain_itail_50nA), 'Color', color_itail_50nA);

set(gca, 'XScale', 'log'); % Set the x-axis as log scale

xlim([1,10e3])
ylim([0.1*min(20*log10(mean_gain)),1.5*max(20*log10(mean_gain))])
%legend('Extrapolated A_{CL} (from measured A_{OL})','Measured A_{OL}','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);


saveas(gcf,[baseName,'multitrials-points-combinedacrossitail','.fig'])
saveas(gcf,[baseName,'multitrials-points-combinedacrossitail','.png'])

%% COMBINED ACROSS ITAILS - extrapolated only
figure
% Define the new colors using the RGB triplets
color_itail_1nA = [0 0.4470 0.7410];   % Dark blue
color_itail_10nA = [0.8500 0.3250 0.0980];  % Dark orange
color_itail_50nA = [0.9290 0.6940 0.1250]; % Dark yellow

% Scatter plots with itail_1nA data (set to blue)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_1nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_1nA); hold on;
% Tracing line for itail_1nA data
h1=plot(mean_freq, 20*log10(mean_gain_cl_itail_1nA), 'Color', color_itail_1nA);

% Scatter plots with itail_10nA data (set to green)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_10nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_10nA);
% Tracing line for itail_10nA data
h2=plot(mean_freq, 20*log10(mean_gain_cl_itail_10nA), 'Color', color_itail_10nA);

% Scatter plots with itail_100nA data (set to red)
scatter(mean_freq, 20*log10(mean_gain_cl_itail_50nA), 100, 's', 'filled', 'MarkerFaceColor', color_itail_50nA);
% Tracing line for itail_100nA data
h3=plot(mean_freq, 20*log10(mean_gain_cl_itail_50nA), 'Color', color_itail_50nA);

set(gca, 'XScale', 'log'); % Set the x-axis as log scale

xlim([1,10e3])
ylim([0,25])
legend([h1 h2 h3],'I_{TAIL}=1nA','I_{TAIL}=10nA','I_{TAIL}=50nA','location','northwest')

% Set labels and title
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
%title(['Gain Plot']);

hold off;
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);


saveas(gcf,[baseName,'multitrials-points-combinedacrossitail-extrapolatedonly','.fig'])
saveas(gcf,[baseName,'multitrials-points-combinedacrossitail-extrapolatedonly','.png'])
%% COMBINE IN ONE FIGURE (subplot)
data1 = ax_itail_1nA.Children
props1 = get(ax_itail_1nA)

data2 = ax_itail_10nA.Children
props1 = get(ax_itail_10nA)


figure
subplot(1,3,1)
copyobj(data1, gca); % Copy saved data to subplot
set(gca, props1);    % Apply saved axes properties

subplot(1,3,2)
copyobj(data2, gca); % Copy saved data to subplot
set(gca, props2);    % Apply saved axes properties


%%
fig1 = openfig('output\ID_sample_I2_covered_for_paper_100mVpp_nolight_itail_1nA_bodeplot__gain_vs_fin_sweep-multitrials-boxplot.fig')
fig2 = openfig('output\ID_sample_I2_covered_for_paper_100mVpp_nolight_itail_10nA_bodeplot__gain_vs_fin_sweep-multitrials-boxplot.fig')
fig3 = openfig('output\ID_sample_I2_covered_for_paper_100mVpp_nolight_itail_50nA_bodeplot__gain_vs_fin_sweep-multitrials-boxplot.fig')

figure
subplot(1,3,1)
ax1 = gca;
copyobj(allchild(get(fig1, 'CurrentAxes')), ax1);

subplot(1,3,2)
ax2 = gca;
copyobj(allchild(get(fig2, 'CurrentAxes')), ax2);

subplot(1,3,3)
ax3 = gca;
copyobj(allchild(get(fig3, 'CurrentAxes')), ax3);
%%
clear all

