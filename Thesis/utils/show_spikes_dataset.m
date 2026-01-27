folderPath = pwd;
runID='intuitivechip_classifier_lcsvsnlcs_teensy_run2'
% name='sensor_1b_pdms_N5b_itail_50nA_squ_indentperiod-trunc-1600-sorted'
% name='N_4b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_180mvpp_250mvofs_onedac'
%name='IMODE_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_9mvpp_250mvofs_onedac_itail_50nA_ith_9p375n_ispk_100n_islew_1n_vil_0p7'
%name='IMODE_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_9mvpp_250mvofs_onedac_itail_50nA_ith_10n_ispk_50u_islew_1n_vil_0p775'
name='IMODE_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_9mvpp_250mvofs_onedac_itail_50nA_ith_9p375n_ispk_100n_islew_1n_vil_0p8' 
name=['VMODE_allchannels_ALLCLASS_auto_corr_cleaned_1x8_18mvpp_250mvofs_onedac_itail_50nA_N5b']
%name='VMODE_allchannels_ALLCLASS_auto_nospkfilter_cleaned_1x8_18mvpp_250mvofs_onedac_itail_50nA_N3b_OFFmerg'
name='IMODE_allchannels_ALLCLASS_auto_corr_cleaned_1x8_180mvpp_250mvofs_onedac_itail_50nA_ith_10n_ispk_50u_islew_1n_vil_1p1'
%name='IMODE_allchannels_ALLCLASS_auto_nospkfilter_cleaned_1x8_180mvpp_250mvofs_onedac_itail_50nA_ith_9p375n_ispk_1p4u_islew_1n_vil_0p9'

tex_cnt=6 %6;  % if <108, choose only a subset of the dataset for easier classification
trialcnt=10 %40  % how many trials per class to use
mr_cnt=8%12*16;


subfigcnt=1;
figure


spike_count_tot=0;
for k=1:tex_cnt
    for trial=1:trialcnt

        filename = [folderPath,'/',runID,'/SENSOR_DATASETS/PROCESSED_SPIKES/',name,'/',num2str(trial+(k-1)*40),'.bs2']
        TD = Read_Ndataset(filename);
        %TD = Read_Ndataset_longTs(filename);
        % check if non-monotically increasing timestamps
        if any(diff(TD.ts<0))
            error('Non-monotonic timestamp!')
        end
        
        spike_count_tot = spike_count_tot + length(TD.ts);

        subplot(tex_cnt,trialcnt,subfigcnt)
        subfigcnt=subfigcnt+1

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
        
        xlim([0 400])
        ylim([0 mr_cnt-1])
        grid on
        xlabel('Time (ms)');
        ylabel('Taxel Number');
        %title([folderPath,'/out/',name,'/',num2str(trial),append,'.bs2']);
        
        TD_on = [];
        TD_off = [];
        
        TD_on.x = TD.x(TD.p==1);
        TD_on.y = TD.y(TD.p==1);
        TD_on.p = TD.p(TD.p==1);
        TD_on.ts = TD.ts(TD.p==1);
        
        TD_off.x = TD.x(TD.p==2);
        TD_off.y = TD.y(TD.p==2);
        TD_off.p = TD.p(TD.p==2);
        TD_off.ts = TD.ts(TD.p==2);
        
        ISI_on = diff(TD_on.ts);
        ISI_off = diff(TD_off.ts);
        
        %micasplot
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
    end
end

spike_count_tot
spike_count_ave = spike_count_tot / (tex_cnt * trialcnt)

sgtitle([name,'-',num2str(spike_count_ave)])
