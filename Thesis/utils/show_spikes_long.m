%% USED TO PLOT SPIKES DIRECTLY GATHERED THROUGH MATLAB (read_serial_spikes.m)
%close all

folderPath = pwd;
name = 'sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod_50ms-2023-08-29_01-31-53_sorted';
name='sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod'
trial=23;
chunks=0;  % 0 if dataset without the -x


if chunks==0
    append = '';
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/SENSOR_DATASETS/PROCESSED_SPIKES/',name,'-long/',num2str(trial),append,'.bs2']
%filename = [folderPath,'/output_datasets\sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod_/',num2str(trial),append,'.bs2']
%TD = Read_Ndataset(filename);
TD = Read_Ndataset_longTs(filename);

figure
hold all;
for spikeCount = 1:length(TD.ts)
    if TD.p(spikeCount)==1  % ON spike
        plot([TD.ts(spikeCount)/1e6 TD.ts(spikeCount)/1e6], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
    else   % OFF spike
       plot([TD.ts(spikeCount)/1e6 TD.ts(spikeCount)/1e6], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
    end 
end

grid on
xlabel('Time (s)');
ylabel('Taxel Number');
title([folderPath,'/out/',name,'/',num2str(trial),append,'.bs2']);

micasplot

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

% ISI_on = diff(TD_on.ts)
% ISI_off = diff(TD_off.ts)
% 
% figure
% % histogram(ISI_on,'BinWidth',5);
% histogram(ISI_on)
% 
% grid on
% xlabel('ISI on spikes (ms)');
% 
% micasplot
% 
% figure
% % histogram(ISI_off,'BinWidth',5);
% histogram(ISI_off)
% 
% grid on
% xlabel('ISI off spikes (ms)');
% size(TD.ts)
% 
% micasplot

