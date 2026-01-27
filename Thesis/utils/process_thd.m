%% READ output of keysight_rigol_thd.py over N iterations for averaging
close all
run_iter = 1;

%% ---- ITAIL=1nA
run_id = 'sample_I2_covered_for_paper_itail_1nA'
baseName = ['output/ID_', run_id, '_thd_vamp_sweep__thd_vs_vamp_sweep-']

vofs=-0.06;

vamp=[]; vout_thd=[]; vin_thd=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csvFileName
    csv_val = csvread(csvFileName,2,0); % skip first line of names & bad first value
    
    vamp(i,:) = csv_val(:,1);
    vout_thd(i,:) = csv_val(:,2);
    vin_thd(i,:) = csv_val(:,3);

end    

vamp=vamp'/100; % attenuation 100x
vout_thd=vout_thd';
vin_thd=vin_thd';

figure
plot(1e3*vamp,vout_thd, 'MarkerFaceColor', 'b'); 
hold on;
scatter(1e3*vamp, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input amplitude (mVpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
hold off;
%--- PLOT as Charge input
Cf=33e-15;
Av=140;
Csensor=172e-15;
chg_tf=(Csensor + Av*Cf); % referred as charge
Cf_max=840e-15;
chg_tf_cf_max=(Csensor + Av*Cf_max); % referred as charge

vcharge=vamp.*chg_tf;
vcharge_cf_max = vamp.*chg_tf_cf_max;
figure
subplot(2,1,1)
plot(1e15*vcharge,vout_thd);
hold on;
scatter(1e15*vcharge, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');
grid on;
%xlabel('Charge-referred input (fC_{pp})');
ylabel('THD (%)');
ylim([0 1.2*max(vout_thd)])
title('Gain setting = 00 (C_{f}=40fF)')

subplot(2,1,2)
plot(1e15*vcharge_cf_max,vout_thd);
hold on;
scatter(1e15*vcharge_cf_max, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Charge-referred input (fC_{pp})');
ylabel('THD (%)');
ylim([0 1.2*max(vout_thd)])
title('Gain setting = 11 (C_{f}=840fF)')

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 400;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-chg','.fig'])
saveas(gcf,[baseName,'multitrials-chg','.png'])



%% ---- ITAIL=10nA
run_id = 'sample_I2_covered_for_paper_itail_10nA'
baseName = ['output\ID_', run_id, '_thd_vamp_sweep__thd_vs_vamp_sweep-']

vofs=-0.02;

vamp=[]; vout_thd=[]; vin_thd=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csvFileName
    csv_val = csvread(csvFileName,2,0); % skip first line of names & bad first value
    
    vamp(i,:) = csv_val(:,1);
    vout_thd(i,:) = csv_val(:,2);
    vin_thd(i,:) = csv_val(:,3);

end    

vamp=vamp'/100; % attenuation 100x
vout_thd=vout_thd';
vin_thd=vin_thd';

figure
plot(1e3*vamp,vout_thd, 'MarkerFaceColor', 'b'); 
hold on;
scatter(1e3*vamp, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input amplitude (mVpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
hold off;
%--- PLOT as Charge input
Cf=33e-15;
Av=140;
Csensor=172e-15;
chg_tf=(Csensor + Av*Cf); % referred as charge

vcharge=vamp.*chg_tf;
figure
plot(1e15*vcharge,vout_thd);
hold on;
scatter(1e15*vcharge, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input charge (fCpp)');
ylabel('THD (%)');

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-chg','.fig'])
saveas(gcf,[baseName,'multitrials-chg','.png'])


%% ---- ITAIL=100nA
run_id = 'sample_I2_covered_for_paper_itail_100nA'
baseName = ['output\ID_', run_id, '_thd_vamp_sweep__thd_vs_vamp_sweep-']

vofs=-0.12;

vamp=[]; vout_thd=[]; vin_thd=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csvFileName
    csv_val = csvread(csvFileName,2,0); % skip first line of names & bad first value
    
    vamp(i,:) = csv_val(:,1);
    vout_thd(i,:) = csv_val(:,2);
    vin_thd(i,:) = csv_val(:,3);

end    

vamp=vamp'/100; % attenuation 100x
vout_thd=vout_thd';
vin_thd=vin_thd';

figure
plot(1e3*vamp,vout_thd, 'MarkerFaceColor', 'b'); 
hold on;
scatter(1e3*vamp, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input amplitude (mVpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
hold off;
%--- PLOT as Charge input
Cf=33e-15;
Av=140;
Csensor=172e-15;
chg_tf=(Csensor + Av*Cf); % referred as charge

vcharge=vamp.*chg_tf;
figure
plot(1e15*vcharge,vout_thd);
hold on;
scatter(1e15*vcharge, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input charge (fCpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-chg','.fig'])
saveas(gcf,[baseName,'multitrials-chg','.png'])

%% ---- ITAIL=50nA
run_id = 'sample_I2_covered_for_paper_itail_50nA'
baseName = ['output\ID_', run_id, '_thd_vamp_sweep__thd_vs_vamp_sweep-']

vofs=-0.12;

vamp=[]; vout_thd=[]; vin_thd=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csvFileName
    csv_val = csvread(csvFileName,2,0); % skip first line of names & bad first value
    
    vamp(i,:) = csv_val(:,1);
    vout_thd(i,:) = csv_val(:,2);
    vin_thd(i,:) = csv_val(:,3);

end    

vamp=vamp'/100; % attenuation 100x
vout_thd=vout_thd';
vin_thd=vin_thd';

figure
plot(1e3*vamp,vout_thd, 'MarkerFaceColor', 'b'); 
hold on;
scatter(1e3*vamp, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input amplitude (mVpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
hold off;
%--- PLOT as Charge input
Cf=33e-15;
Av=140;
Csensor=172e-15;
chg_tf=(Csensor + Av*Cf); % referred as charge

vcharge=vamp.*chg_tf;
figure
plot(1e15*vcharge,vout_thd);
hold on;
scatter(1e15*vcharge, vout_thd, 100, 's', 'filled', 'MarkerFaceColor', 'b');

% Set labels and title
xlabel('Input charge (fCpp)');
ylabel('THD (%)');

% Add legend
%legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southwest');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-chg','.fig'])
saveas(gcf,[baseName,'multitrials-chg','.png'])


clear all
