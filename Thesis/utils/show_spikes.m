folderPath = pwd;
% name = 'sensorized_test_1b_pdms_N4p5b_squaretip_indent_sorted';
% name= 'sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod_noofs'
%name='isi_10hz_100mvpp_itail_1nA_sinewave_new'
%name = 'N_3b_allchannels_cl1_auto_nospkfilter_12x16_180mvpp_250mvofs_sorted/oneperiod';
name='IMODE_allchannels_ALLCLASS_auto_corr_cleaned_1x8_180mvpp_250mvofs_onedac_itail_50nA_ith_10n_ispk_50u_islew_1n_vil_0p8'
trial=210;
chunks=0;  % 0 if dataset without the -x
%name = 'latiti'
%name = 'N_3b_allchannels_cl1_auto_nospkfilter_12x16_180mvpp_250mvofs_sorted/oneperiod';

if chunks==0
    append = ''
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/SENSOR_DATASETS/PROCESSED_SPIKES/',name,'/',num2str(trial),append,'.bs2']
TD = Read_Ndataset(filename);
%TD = Read_Ndataset_longTs(filename);

figure
hold all;
for spikeCount = 1:length(TD.ts)
    if TD.p(spikeCount)==1  % ON spike
        plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
    else   % OFF spike
       plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
    end 
end

grid on
xlabel('Time (ms)');
ylabel('Taxel Number');
title([folderPath,'/out/',name,'/',num2str(trial),append,'.bs2']);

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

ISI_on = diff(TD_on.ts)
ISI_off = diff(TD_off.ts)

micasplot
% 
% figure
% histogram(ISI_on,'BinWidth',5)
% 
% grid on
% xlabel('ISI on spikes (ms)');
% 
% figure
% histogram(ISI_off,'BinWidth',5)
% 
% grid on
% xlabel('ISI off spikes (ms)');
