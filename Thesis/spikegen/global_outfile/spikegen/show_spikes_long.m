%% Author: Mark Daniel Alea
% mark.alea@kuleuven.be
%
% reads .bs2 files and display as spatiotemporal spikes
% _long suffix means that the .bs2 file stores spike times as a 31bit
% integer (microsecond resolution) -> especially useful for long spike
% recordings

folderPath = pwd;
name = 'N_4b_allchannels_cl1_auto_nospkfilter_12x16_180mvpp_250mvofs_onedac';
name = 'N_4b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_180mvpp_250mvofs_onedac'
%name = 'ALL_N_lcadc-0-uV-v2-0-uV-0-uV-0-pct-4-mr_cnt-192-multitrial-Kylberg_filt_6_scan_actualdimscale_chip_lchpf--1to40-1to40_sorted'
%name = 'sensorized_test_1a_0_n5b_squaretip_xyz_indent_pdmsreal_short_trial3_sorted'
trial=1;
chunks=0;  % 0 if dataset without the -x
%name = 'latiti'
%name = 'N_3b_allchannels_1p80mvpp_250mvofs_cl1_tr1';

if chunks==0
    append = '';
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/global_outfile/spikegen/',name,'/',num2str(trial),append,'.bs2']
filename = ['/users/micas/malea/Documents/ArduinoRead_2/ArduinoRead/out/',name,'/',num2str(trial),append,'.bs2']
TD = Read_Ndataset(filename);
%TD = Read_Ndataset_longTs(filename);

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
micasplot
