function [hem] = APDSM_sigma_Quantization(current_loc,tsel,trialsel,N,LSB,LSB_th,synfilt_en,tau_syn,tau_refr,stim,max_maxSR,MR_type,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,noiseVar,th_mismatch,tau_refr_mismatch,scandir,fileprint)

%clear all
close all

tic

tex_rs = stim.tex_rs;
tex_rsdiff = stim.tex_rsdiff;
tex_name = stim.tex_name;
L = stim.L;
t = stim.t;
ti = stim.ti;
tlength = stim.tlength;
fs_stim = stim.fs_stim;
stim_os = stim.stim_os;
DT=stim.DT; % needed for tau computation

%DT=1/(fs_stim*stim_os);
%DT=(scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
disp(' '); disp(' '); 
disp(['APDSM - Processing texture: ',tex_name{tsel},'  trial # ',num2str(trialsel)]);
disp(['DT=',num2str(DT)]);

% sigbw=100;
% osr=1;
% fs=2*sigbw*osr;

%sizeof sampled matrix = 20secs * fs
%sizeof resampled stimulus = 20secs * 200 *stim_os

%% APDSM quantisation
%N = 2;
%LSB = 1/(2^N);
%tex_rs_chunk = tex_rs{tsel}(n,:); % get 1 second chunk
if MR_type == 'FA'
    tex_rs_chunk = tex_rsdiff{tsel,trialsel}; % get whole chunk
elseif MR_type == 'SA'
    tex_rs_chunk = tex_rs{tsel,trialsel}; % get whole chunk   
else
    tex_rs_chunk = 0; % get whole chunk
end

if crf==1
    cnt=hem_cnt;
else
    cnt=mr_cnt;
end


% %% MODEL low-freq variability trial-to-trial (assumes each heminode has independent trial-to-trial variation)
% trialNoise = randn(cnt,length(t));
% %lpFilter = designfilt('lowpassiir','FilterOrder',10,'HalfPowerFrequency1',100,'SampleRate',500);
% lpFilter = designfilt('lowpassiir', 'FilterOrder', 10, ...
%                       'HalfPowerFrequency', 0.5, 'SampleRate', 500);
% for h=1:cnt                 
%     trialNoise_filt(h,:) = filter(lpFilter,trialNoise(h,:));
%     trialNoise_filt_scaled(h,:) = trialNoise_filt(h,:).*(noiseVar(1)/(max(trialNoise_filt(h,:)) - min(trialNoise_filt(h,:))));
% end

%% MODEL low-freq variability trial-to-trial 
if noiseVar(1) > 0    % only do expensive processing if needed
    trialNoise = randn(1,stim.tlength_origrs);
    %lpFilter = designfilt('lowpassiir','FilterOrder',10,'HalfPowerFrequency1',100,'SampleRate',500);
    lpFilter = designfilt('lowpassiir', 'FilterOrder', 10, ...
                          'HalfPowerFrequency', 0.5, 'SampleRate', 500);

    trialNoise_filt = filter(lpFilter,trialNoise);
    trialNoise_filt_scaled_ = trialNoise_filt.*(noiseVar(1)/(max(trialNoise_filt) - min(trialNoise_filt)));
    temp = reshape(trialNoise_filt_scaled_ ,length(trialNoise_filt_scaled_)/stim.time_chunks,stim.time_chunks); % divide texture sample to groups of tdelay seconds
    temp = temp';
    for tt=1:stim.tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
        %disp(['--->Creating time chunk # ',num2str(t)]);  
        temp_t = [];
        for i=1:round(stim.tlength_des/stim.tdelay)
            temp_t = [temp_t,temp(i+tt-1,:)];
        end
        trialNoise_stitched{tt} = temp_t;   % different heminodes get different 1s section (DELAY IN OVERLAP WITH CRF)
    end
% 
%     td = stim.tdelaycnt; 
%     for h=0:mr_cnt/8 - 1
%         %disp(['--->Filling heminodes # ',num2str(8*h+1),'-',num2str(8*h+8)]);
%         trialNoise_filt_scaled(8*h+1,:) = trialNoise_stitched{td};
%         trialNoise_filt_scaled(8*h+2,:) = trialNoise_stitched{td-1};            
%         trialNoise_filt_scaled(8*h+3,:) = trialNoise_stitched{td-2};
%         trialNoise_filt_scaled(8*h+4,:) = trialNoise_stitched{td-3};         
%         trialNoise_filt_scaled(8*h+5,:) = trialNoise_stitched{td-4};
%         trialNoise_filt_scaled(8*h+6,:) = trialNoise_stitched{td-5};            
%         trialNoise_filt_scaled(8*h+7,:) = trialNoise_stitched{td-6};
%         trialNoise_filt_scaled(8*h+8,:) = trialNoise_stitched{td-7};   
%         td = td-1;
%     end
    [trialNoise_filt_scaled] = fill_array(trialNoise_stitched,1,1,stim.tdelaycnt,mr_cnt,scandir); 
     


else
    trialNoise = randn(1,stim.tlength_origrs);    
    trialNoise_filt_scaled = zeros(mr_cnt,length(t));
end

%% MODEL sensor cross-talk
sensorxtalk = sqrt(noiseVar(3))*randn(mr_cnt,length(t));
%% MODEL threshold mismatch (DC)
%th_mismatch = sqrt(noiseVar(4))*randn(1,mr_cnt);
%% threshold noise (AC) + mismatch (DC) modeled inside quantization routine
th_noise = sqrt(noiseVar(2))*randn(cnt,length(t));   % generate Gaussian threshold noise (unique for every heminode, every trial, etc)

ff43=figure(43);
subplot(4,1,1)
plot(t,trialNoise(1:length(t)),'linewidth',2)
title('Trial variability (unfiltered)')
subplot(4,1,2)
plot(t,trialNoise_filt_scaled(1,:),'linewidth',2); hold on;
%plot(t,trialNoise_filt2(1,:),'linewidth',2);
title('Trial variability (filtered <0.5Hz + scaled)')
subplot(4,1,3)
plot(t,sensorxtalk(1,:),'linewidth',2)
title('Sensor xtalk')
subplot(4,1,4)
plot(t,th_noise(1,:) + th_mismatch(1),'linewidth',2)
title('Threshold noise + mismatch')

ff44=figure(44);
subplot(4,1,1)
plot(t,trialNoise(1:length(t)),'linewidth',2)
title('Trial variability (unfiltered)')
subplot(4,1,2)
plot(t,trialNoise_filt_scaled(cnt,:),'linewidth',2); hold on;
%plot(t,trialNoise_filt2(1,:),'linewidth',2);
title('Trial variability (filtered <0.5Hz + scaled)')
subplot(4,1,3)
plot(t,sensorxtalk(cnt,:),'linewidth',2)
title('Sensor xtalk')
subplot(4,1,4)
plot(t,th_noise(cnt,:) + th_mismatch(cnt),'linewidth',2)
title('Threshold noise + mismatch')

tex_rs_chunk = tex_rs_chunk + trialNoise_filt_scaled + sensorxtalk;

[u,q,q_node,delta,vth,g_syn,y,spktimes,spktimes_node] = APDSM_sigma(tex_rs_chunk,t,max_maxSR,N,DT,LSB,LSB_th,synfilt_en,tau_syn,tau_refr,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,th_noise,th_mismatch,tau_refr_mismatch);
size(y)
if wta == 1
    hem.spktimes = spktimes_node;
else
    hem.spktimes = spktimes;
end

%% figures
% ff1=figure(1);
% subplot(4,1,1)
% plot(t,tex_rs_chunk,'linewidth',2);
% title('Original signal');
% xlabel('Time');
% ylabel('Amplitude');
% subplot(4,1,2);
% stem(t,q);
% title(['Q signal (N = ',num2str(N),'; LSB = ',num2str(LSB),' mV)']);
% xlabel('Time');
% ylabel('Amplitude');
% subplot(4,1,3);
% stem(t,dir);
% title(['Dir signal (N = ',num2str(N),'; LSB = ',num2str(LSB),' mV)']);
% xlabel('Time');
% ylabel('Amplitude');

%% reconstruction algorithm
% w = 2*pi*sigbw/osr;
% z = 0;
% for k = 1:length(t)
%     %disp(['reconstructing index',num2str(k),' - ',num2str(100*k/length(t)),'%']);
%     z = z + LSB.*dir(k);
% end 
% c = max(tex_rs_chunk)./max(z); %scaling is done as a consequence of oversampling
% z = z.*c;

%% figures
% subplot(4,1,4)
% plot(t,u,'linewidth',2);
% title('Reconstructed signal');
% xlabel('Time');
% ylabel('Amplitude');

ff2=figure(12);
clf(ff2);
subplot(3,1,1)
plot(t,y(1,:),'linewidth',2)
hold on
plot(t,u(1,:),'linewidth',2);
plot(t,vth(1,:),'linewidth',2);
%plot(t,y_sigma,'linewidth',2);
title('APDSM Quantization');
legend('Input signal','u','vth');
%ylim([0 1]);
xlim([0 max(t)]);
subplot(3,1,2)
plot(t,delta(1,:),'linewidth',2);
line([0 4], [LSB_th/2+th_mismatch(1) LSB_th/2+th_mismatch(1)],'LineStyle','--','Color','r');
legend('delta','vth=LSB_th/2');
xlim([0 max(t)]);
ylim([-LSB_th LSB_th]);
subplot(3,1,3)
stem(t,q(1,:),'linewidth',2);
legend('q');
xlim([0 max(t)]);

if cnt > 1
    ff22=figure(22);
    clf(ff22);
    subplot(3,1,1)
    plot(t,y(2,:),'linewidth',2)
    hold on
    plot(t,u(2,:),'linewidth',2);
    plot(t,vth(2,:),'linewidth',2);
    %plot(t,y_sigma,'linewidth',2);
    title('APDSM Quantization');
    legend('Input signal','u','vth');
    %ylim([0 1]);
    xlim([0 max(t)]);
    subplot(3,1,2)
    plot(t,delta(2,:),'linewidth',2);
    line([0 4], [LSB_th/2+th_mismatch(2) LSB_th/2+th_mismatch(2)],'LineStyle','--','Color','r');
    legend('delta','vth = LSB_th/2');
    xlim([0 max(t)]);
    ylim([-LSB_th LSB_th]);
    subplot(3,1,3)
    stem(t,q(2,:),'linewidth',2);
    legend('q');
    xlim([0 max(t)]);
end

if crf==1
    cnt=hem_cnt;
else
    if mr_cnt > 64
        %cnt=64;  % fit as many heminodes
        cnt=mr_cnt;
    else
        cnt=mr_cnt;
    end
end

ff23=figure(23);
clf(ff23);
set(gcf,'position',[10,10,1800,1200])
for h=1:cnt
    %subplot(cnt+1,1,h)
    subplot(sqrt(mr_cnt)+1,sqrt(mr_cnt),h);
    stem(t,q(h,:),'linewidth',2);  
    xlim([0 max(t)]);
    title(['Heminode #',num2str(h)])
end
%subplot(cnt+1,1,cnt+1)
subplot(sqrt(mr_cnt)+1,sqrt(mr_cnt),cnt+1)
stem(t,q_node(:),'linewidth',2); 
title('Spikes from node of Ranvier (WTA)')


ff24=figure(24);
clf(ff24);
set(gcf,'position',[10,10,1800,1200])
for h=1:cnt
    %subplot(cnt+1,1,h)
    subplot(sqrt(mr_cnt),sqrt(mr_cnt),h);
    plot(t,y(h,:),'linewidth',2);  
    xlim([0 max(t)]);
    title(['Heminode #',num2str(h)])
end

for h=1:cnt
    SR{h} = gradient(y(h,:))/(t(2)-t(1));
    maxSR(h) = max(abs(SR{h}));   
    tlooprequired(h) = LSB(h)/maxSR(h);
    disp(['Required tloop to avoid slope overload: <',num2str(tlooprequired(h))]);
end
%hem.tlooprequired = tlooprequired;
%hem.maxSR = maxSR;

ff3=figure(3);
plot(t,SR{1},'linewidth',2)
title('Slew rate of 1st heminode')


%% mse computation
for h=1:cnt
    error(h) = immse(u(h,:),y(h,:));
    disp(['MSE for hem #',num2str(h),' = ',num2str(error(h))]);
    %Psig = rms(tex_rs{tsel}(n,:))^2
    Psig(h) = rms(y(h,:))^2;
    disp(['Psig for hem #',num2str(h),' = ',num2str(Psig(h))]);
    SER(h) = 10*log10(Psig(h)/error(h));
    disp(['SER (dB) for hem #',num2str(h),' = ',num2str(SER(h))]);
end
hem.SER = SER;

%% get spike information
for h=1:cnt
    %spikecount(h) = length(find(q(h,:)==1));
    spikecount(h) = numel(spktimes{h});
    disp(['for heminode # ',num2str(h),': spikerate = ',num2str(spikecount(h)/tlength),' Hz']); 
end
hem.spikecount = spikecount;

%% display number of samples for all heminodes
for h=1:cnt
    %samplecnt(h) = length(find(q(h,:)==1));
    samplecnt(h) = numel(spktimes{h});
    disp(['APDSM   N=',num2str(N),'  LSB= ',num2str(LSB(h)),'  TH = ',num2str(LSB_th+th_mismatch(h)),'  FS tau_refr=',num2str(tau_refr+tau_refr_mismatch(h)),'   texture:',tex_name{tsel},'  samples for hem #',num2str(h),' = ',num2str(samplecnt(h))]);
end
hem.samplecnt = samplecnt;

if fileprint==1
    % fname1 = (['figures/APDSMtexture-',tex_name{tsel},'-','s_fig1']);
     fname2 = ([current_loc,'figures/spikegen/APDSMtextureAR-',tex_name{tsel},'- N=',num2str(N),'s_fig2']);
    % fname3 = (['figures/APDSMtexture-',tex_name{tsel},'-','s_fig3']);
    % fname4 = (['figures/APDSMtexture-',tex_name{tsel},'-','s_fig4']);
    % print(ff1,'-dpng','-r600',fname1);
     print(ff2,'-dpng','-r600',fname2);
    % print(ff3,'-dpng','-r600',fname3);
    % print(ff4,'-dpng','-r600',fname4);
end
toc

end
