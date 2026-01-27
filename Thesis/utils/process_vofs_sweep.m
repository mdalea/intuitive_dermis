%% READ output of keysight_rigol_ofs_sweep.py over N iterations for averaging
run_iter = 3;

%% ---- ITAIL=1nA
run_id = 'sample_I2_covered_for_paper_itail_1nA'
%run_id = 'sample_I2_covered_for_paper_100mVpp_nolight_itail_1nA'
%run_id = 'sample_I2_itail_100nA'
baseName = ['output\ID_', run_id, '_ofs_sweep__gain_vs_vofsin_sweep-']

% used for TSENSE
vamp=0.05;
fin=10;
vofs_arr = [-200e-3:20e-3:200e-3];

vofs=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    vofs(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

vofs=vofs'; vofs=-1*vofs/100; % 100x attenuation and inversion
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

figure
plot(vofs,gain)

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southeast');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
%---------
figure
mean_vofs = mean(vofs,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get optimum vofs where gain is at max
vofs_opt = mean_vofs(find(mean_gain==max(mean_gain)));

% Plot the means as a solid line
plot(mean_vofs, mean_gain, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_vofs, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','northwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = mean_vofs(round(length(mean_vofs) / 2));
y_pos = min(mean_gain);
text(x_pos, y_pos, sprintf('vofs_{opt} = %.6f V', vofs_opt), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])
% Plot temperature over iterations and sweep values
baseName = ['output\ID_', run_id, '_ofs_sweep__tsense_ave_fin_', num2str(fin),'_vpp_', num2str(vamp), '_vofs_']


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(vofs_arr)
        if vofs_arr(j) == 0            
            csvFileName = [baseName,'0.0','-',num2str(i-1),'.csv']; % that pesky case where num2str() gets confused with python's str()
        else
            csvFileName = [baseName,num2str(vofs_arr(j)),'-',num2str(i-1),'.csv'];
        end    
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
vofs_arr=vofs_arr'; vofs_arr=-1*vofs_arr/100 %100x attenuation and inversion
figure
hold on;
plot(vofs_arr,tsense_arr)
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
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
%run_id = 'sample_I2_itail_100nA'
baseName = ['output\ID_', run_id, '_ofs_sweep__gain_vs_vofsin_sweep-']

% used for TSENSE
vamp=0.05;
fin=10;
vofs_arr = [-200e-3:20e-3:200e-3];

vofs=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    vofs(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

vofs=vofs'; vofs=-1*vofs/100; % 100x attenuation and inversion
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

figure
plot(vofs,gain)

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southeast');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
%---------
figure
mean_vofs = mean(vofs,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get optimum vofs where gain is at max
vofs_opt = mean_vofs(find(mean_gain==max(mean_gain)));

% Plot the means as a solid line
plot(mean_vofs, mean_gain, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_vofs, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','northwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = mean_vofs(round(length(mean_vofs) / 2));
y_pos = min(mean_gain);
text(x_pos, y_pos, sprintf('vofs_{opt} = %.6f V', vofs_opt), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])
% Plot temperature over iterations and sweep values
baseName = ['output\ID_', run_id, '_ofs_sweep__tsense_ave_fin_', num2str(fin),'_vpp_', num2str(vamp), '_vofs_']


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(vofs_arr)
        if vofs_arr(j) == 0            
            csvFileName = [baseName,'0.0','-',num2str(i-1),'.csv']; % that pesky case where num2str() gets confused with python's str()
        else
            csvFileName = [baseName,num2str(vofs_arr(j)),'-',num2str(i-1),'.csv'];
        end    
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
vofs_arr=vofs_arr'; vofs_arr=-1*vofs_arr/100 %100x attenuation and inversion
figure
hold on;
plot(vofs_arr,tsense_arr)
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
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
%run_id = 'sample_I2_itail_100nA'
baseName = ['output\ID_', run_id, '_ofs_sweep__gain_vs_vofsin_sweep-']

% used for TSENSE
vamp=0.05;
fin=10;
vofs_arr = [-200e-3:20e-3:200e-3];

vofs=[]; vout_pp=[]; vin_pp=[]; gain=[]; gain_db=[];
for i = 1:run_iter
    csvFileName = [baseName,num2str(i-1),'.csv' ];
    csv_val = csvread(csvFileName,1,0); % skip first line of names
    
    vofs(i,:) = csv_val(:,1);
    vout_pp(i,:) = csv_val(:,2);
    vin_pp(i,:) = csv_val(:,3);
    gain(i,:) = csv_val(:,4);
    gain_db(i,:) = csv_val(:,5);
    
end    

vofs=vofs'; vofs=-1*vofs/100; % 100x attenuation and inversion
vout_pp=vout_pp';
vin_pp=vin_pp';
gain=gain';
gain_db=gain_db';

figure
plot(vofs,gain)

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Gain Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5', 'location', 'southeast');
micasplot

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials','.fig'])
saveas(gcf,[baseName,'multitrials','.png'])
%---------
figure
mean_vofs = mean(vofs,2);  % get mean across dim=2 just to have it 1 column; 
mean_gain = mean(gain,2); % get average across iterations 2nd dim
std_gain = std(gain,0,2); % get sigma

% get optimum vofs where gain is at max
vofs_opt = mean_vofs(find(mean_gain==max(mean_gain)));

% Plot the means as a solid line
plot(mean_vofs, mean_gain, 'b', 'LineWidth', 2);
hold on;

% Plot the sigmas as error bars (standard deviation)
errorbar(mean_vofs, mean_gain, std_gain, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
ylabel('Gain (V/V)');
title(['Mean and Standard Deviation Across ',num2str(run_iter),' Iterations']);

% Add legend
legend('Mean', 'Standard Deviation','location','northwest');

% Add text annotation for vofs_opt at the bottom center
x_pos = mean_vofs(round(length(mean_vofs) / 2));
y_pos = min(mean_gain);
text(x_pos, y_pos, sprintf('vofs_{opt} = %.6f V', vofs_opt), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

micasplot
hold off;

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,[baseName,'multitrials-boxplot','.fig'])
saveas(gcf,[baseName,'multitrials-boxplot','.png'])
% Plot temperature over iterations and sweep values
baseName = ['output\ID_', run_id, '_ofs_sweep__tsense_ave_fin_', num2str(fin),'_vpp_', num2str(vamp), '_vofs_']


tsense_arr=[]; tsense_std_arr=[];
for i = 1:run_iter
    for j=1:numel(vofs_arr)
        if vofs_arr(j) == 0            
            csvFileName = [baseName,'0.0','-',num2str(i-1),'.csv']; % that pesky case where num2str() gets confused with python's str()
        else
            csvFileName = [baseName,num2str(vofs_arr(j)),'-',num2str(i-1),'.csv'];
        end    
        csv_val = csvread(csvFileName); % skip first line of names

        tsense_arr(i,j) = csv_val(1); % get mean
        tsense_std_arr(i,j) = csv_val(2); % get std

    end
end    

tsense_arr=tsense_arr';
tsense_std_arr=tsense_std_arr';
vofs_arr=vofs_arr'; vofs_arr=-1*vofs_arr/100 %100x attenuation and inversion
figure
hold on;
plot(vofs_arr,tsense_arr)
% Plot the sigmas as error bars (standard deviation)
%errorbar([vofs_arr vofs_arr vofs_arr vofs_arr vofs_arr], tsense_arr, tsense_std_arr, 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

% Set labels and title
xlabel('Offset Voltage (V)');
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