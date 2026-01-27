%% code to extract useful info from the indent pulse dataset
clear all
close all

folderPath = pwd;

name='sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod'

spkcount = []; tmax=[]; spkfreq = [];
for trial=1:200

    chunks=0;  % 0 if dataset without the -x
    
    
    if chunks==0
        append = ''
    else
        append = ['-',num2str(chunks)];
    end
    
    filename = [folderPath,'/SENSOR_DATASETS/PROCESSED_SPIKES/',name,'/',num2str(trial),append,'.bs2']
    TD = Read_Ndataset(filename);
    %TD = Read_Ndataset_longTs(filename);
    
    tmax(trial) = max(TD.ts)
    spkcount(trial) = length(TD.ts)
    spkfreq(trial) = length(TD.ts) / (10e-6*max(TD.ts))

end

figure(50)
subplot(3,1,1)
plot(1:200, tmax);
ylabel('Max spike timestamp (us)')
subplot(3,1,2)
plot(1:200, spkcount);
ylabel('Spike counts')
subplot(3,1,3)
plot(1:200, spkfreq);
ylabel('Spike frequency (Hz)')

figure(51)
plot(1:200, spkfreq);
xlabel('Dataset index')
ylabel(['Spike',newline, 'frequency',newline, '(Hz)'])
[max_spkfreq,index] = max(spkfreq)

micasplot
hold off;
%% bin into windows
trial=161%200;
chunks=0;
if chunks==0
    append = ''
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/SENSOR_DATASETS/PROCESSED_SPIKES/',name,'/',num2str(trial),append,'.bs2']
TD = Read_Ndataset(filename);

% optional filter
ts_min = 1000000; ts_max = 3000000;
TD.x = TD.x(TD.ts > ts_min & TD.ts < ts_max)
TD.y = TD.y(TD.ts > ts_min & TD.ts < ts_max)
TD.p = TD.p(TD.ts > ts_min & TD.ts < ts_max)
TD.ts = TD.ts(TD.ts > ts_min & TD.ts < ts_max)

hold off;

figure(52)
subplot(3,1,1)
hold all;
for spikeCount = 1:length(TD.ts)
    spikeCount
    if TD.p(spikeCount)==1  % ON spike
        plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
    else   % OFF spike
       plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
    end 
end

grid on
%xlabel('Time (ms)');
ylabel(['Taxel',newline,'Number']);
xlim([ts_min/1e3 ts_max/1e3])
micasplot

% BINS
TD_on = [];
TD_off = [];

TD_on.x = TD.x(TD.p==1)
TD_on.y = TD.y(TD.p==1)
TD_on.p = TD.p(TD.p==1)
TD_on.ts = TD.ts(TD.p==1)

TD_off.x = TD.x(TD.p==2)
TD_off.y = TD.y(TD.p==2)
TD_off.p = TD.p(TD.p==2)
TD_off.ts = TD.ts(TD.p==2)

% Define the time window size in microseconds (1 millisecond = 1000 microseconds)
window_size = 10e3; %us

% Calculate the number of time windows needed
num_windows = ceil(max(TD.ts) / window_size);

% Initialize an array to store the counts in each time window
window_counts_on = zeros(1, num_windows);
window_counts_off = zeros(1, num_windows);

% Bin the timestamps into time windows and count the number of timestamps in each window
for i = 1:length(TD_on.ts)
    window_index_on = ceil(TD_on.ts(i) / window_size);
    window_counts_on(window_index_on) = window_counts_on(window_index_on) + 1;
end

% Bin the timestamps into time windows and count the number of timestamps in each window
for i = 1:length(TD_off.ts)
    window_index_off = ceil(TD_off.ts(i) / window_size);
    window_counts_off(window_index_off) = window_counts_off(window_index_off) + 1;
end

% Define the time values for the x-axis (in milliseconds)
time_windows = (1:num_windows) * window_size / 1000;

hold off;
subplot(3,1,2)
% Plot the counts vs. time windows
bar(time_windows, window_counts_on);
%xlabel('Time Windows (ms)');
ylabel('Count');
%title(['ON Spikes - Bin width: ',num2str(window_size*1e-3),'ms']);
t1=text(1100,50,['ON Spikes',newline,'Bin width: ',newline,num2str(window_size*1e-3),'ms'], 'Color', [0.6350 0.0780 0.1840]);
xlim([ts_min/1e3 ts_max/1e3])
micasplot

subplot(3,1,3)
% Plot the counts vs. time windows
bar(time_windows, window_counts_off);
xlabel('Time Windows (ms)');
ylabel('Count');
%title(['ON Spikes - Bin width: ',num2str(window_size*1e-3),'ms']);
t2=text(1100,10,['OFF Spikes',newline,'Bin width: ',newline,num2str(window_size*1e-3),'ms'], 'Color', [0.6350 0.0780 0.1840]);
xlim([ts_min/1e3 ts_max/1e3])
micasplot

set(t1, 'FontSize',10);
set(t2, 'FontSize',10);

% Set the width and height of the figure (in pixels)
width = 800;    % Specify your desired width here
height = 400;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,['SNN_SVM_plots\spike_frames_bin_',name,'-',num2str(trial),'.fig'])
saveas(gcf,['SNN_SVM_plots\spike_frames_bin_',name,'-',num2str(trial),'.png'])
%% heatmap

frame_data_on=zeros(12,16); frame_data_off=zeros(12,16);

figure(53)
subplot(3,1,1)
hold all;
for spikeCount = 1:length(TD.ts)
    spikeCount
    if TD.p(spikeCount)==1  % ON spike
        plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
        %if TD.ts(spikeCount) > 4260000 & TD.ts(spikeCount) < 4280000
            frame_data_on(TD.y(spikeCount)+1,TD.x(spikeCount)+1) = frame_data_on(TD.y(spikeCount)+1,TD.x(spikeCount)+1) + 1 % accumulate spikes per address
        %end
    else   % OFF spike
       plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
       %if TD.ts(spikeCount) > 4260000 & TD.ts(spikeCount) < 4280000
        frame_data_off(TD.y(spikeCount)+1,TD.x(spikeCount)+1) = frame_data_off(TD.y(spikeCount)+1,TD.x(spikeCount)+1) + 1 % accumulate spikes per address
       %end
    end 
end
hold off;

grid on
xlabel('Time (ms)');
ylabel(['Taxel',newline,'Number']);
xlim([0 max(TD.ts)/1e3])
micasplot

subplot(3,1,2)
heatmap(frame_data_on);
colormap(jet(512))
title('ON spike count')
%micasplot

subplot(3,1,3)
heatmap(frame_data_off);
colormap(jet(512))   
title('OFF spike count')
%micasplot

% Set the width and height of the figure (in pixels)
width = 1200;    % Specify your desired width here
height = 600;   % Specify your desired height here
set(gcf, 'Position', [100, 100, width, height]);
saveas(gcf,['SNN_SVM_plots\spike_frames_',name,'-',num2str(trial),'.fig'])
saveas(gcf,['SNN_SVM_plots\spike_frames_',name,'-',num2str(trial),'.png'])