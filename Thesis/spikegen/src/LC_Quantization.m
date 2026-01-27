function [hem] = LC_Quantization(current_loc,tsel,trialsel,N,stim,max_maxSR,MR_type,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,noiseVar,th_mismatch_h,th_mismatch_l,scandir,fileprint)


%clc
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

DT=1/(fs_stim*stim_os);
disp(' '); disp(' '); 
disp(['LC-ADC - Processing texture: ',tex_name{tsel},'  trial # ',num2str(trialsel)]);
disp(['DT=',num2str(DT)]);

sigbw=100;
osr=1;
fs=2*sigbw*osr;


%sizeof sampled matrix = 20secs * fs
%sizeof resampled stimulus = 20secs * 200 *stim_os

%% LC quantisation
%N = 3;
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
    %cnt=hem_cnt;
    cnt=mr_cnt; % preserve original sensor size (bec. here yo uadd all the noise, etc)try for now
else
    cnt=mr_cnt;
end


%% MODEL low-freq variability trial-to-trial (assumes each heminode has independent trial-to-trial variation)
if noiseVar(1) > 0    % only do expensive processing if needed
    trialNoise = randn(1,stim.tlength_origrs);
    %lpFilter = designfilt('lowpassiir','FilterOrder',10,'HalfPowerFrequency1',100,'SampleRate',500);
    lpFilter = designfilt('lowpassiir', 'FilterOrder', 10, ...
                          'HalfPowerFrequency', 0.01, 'SampleRate', 500);

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
%    [trialNoise_filt_scaled] = fill_array(trialNoise_stitched,1,stim.tdelaycnt,mr_cnt); 
    [trialNoise_filt_scaled] = fill_array(trialNoise_stitched,1,1,stim.tdelaycnt,mr_cnt,scandir);
   %commented out!!!!!
    %[trialNoise_filt_scaled] = fill_array_12x16(trialNoise_stitched,1,1,stim.tdelaycnt,mr_cnt,scandir);  


else
    trialNoise = randn(1,stim.tlength_origrs);    
    trialNoise_filt_scaled = zeros(mr_cnt,length(t));
end

%% MODEL AFE white noise
sensorxtalk = sqrt(noiseVar(3))*randn(mr_cnt,length(t));
%% MODEL threshold mismatch (DC)
%th_mismatch = sqrt(noiseVar(4))*randn(1,hem_cnt);
%% threshold noise (AC) + mismatch (DC) modeled inside quantization routine 
th_noise_h = sqrt(noiseVar(2))*randn(cnt,length(t));   % generate Gaussian threshold noise (unique for every heminode, every trial, etc)
th_noise_l = sqrt(noiseVar(2))*randn(cnt,length(t));   % generate Gaussian threshold noise (unique for every heminode, every trial, etc)
%% MODEL AFE flicker noise
length(t)
totalSize=length(t);
df=1;
f=0:1/df:totalSize-1/df;
for titi=1:cnt
    S=1./sqrt(f); 
    S=S.*exp(j*2*pi*randn(size(f))); 
    %S=(8.75/7)*1e3.*S; % normalize
    S=sqrt(noiseVar(5)).*S; % multiplier from 1e-6 V^2 at 1Hz
    S=sqrt(10).*S;
    S(1)=0;
    sensorxtalk_flicker(titi,:) = real(ifft(S));
end    

samplingRate = fs; % Sampling rate (5 MHz in this example)
% windowTime = 1;
% windowSize = samplingRate * windowTime;
% totalSize=windowSize
% overlap = windowSize / 2;
windowSize=length(t);
overlap=0;
[psdEstimate, freq] = pwelch(sensorxtalk_flicker(1,:), windowSize, overlap, windowSize, samplingRate, 'psd', 'onesided'); % psd not power matches better with sim

% Plot the power spectral density
figure;
loglog(freq, psdEstimate); 
grid on;

ff43=figure(43);
subplot(6,1,1)
plot(t,trialNoise(1:length(t)),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Trial variability (unfiltered)')
subplot(6,1,2)
plot(t,trialNoise_filt_scaled(1,:),'linewidth',2); hold on;
%plot(t,trialNoise_filt2(1,:),'linewidth',2);
title('Trial variability (filtered <0.5Hz + scaled)')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,3)
plot(t,sensorxtalk(1,:),'linewidth',2)
title('AFE noise')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,4)
plot(t,sensorxtalk_flicker(1,:),'linewidth',2)
title('AFE noise - flicker')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,5)
plot(t,th_noise_h(1,:) + th_mismatch_h(1),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (upper) + mismatch')
subplot(6,1,6)
plot(t,th_noise_l(1,:) + th_mismatch_l(1),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (lower) + mismatch')

ff44=figure(44);
subplot(6,1,1)
plot(t,trialNoise(1:length(t)),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Trial variability (unfiltered)')
subplot(6,1,2)
plot(t,trialNoise_filt_scaled(cnt,:),'linewidth',2); hold on;
%plot(t,trialNoise_filt2(1,:),'linewidth',2);
title('Trial variability (filtered <0.5Hz + scaled)')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,3)
plot(t,sensorxtalk(cnt,:),'linewidth',2)
title('AFE noise')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,4)
plot(t,sensorxtalk_flicker(cnt,:),'linewidth',2)
title('AFE noise - flicker')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(6,1,5)
plot(t,th_noise_h(cnt,:) + th_mismatch_h(cnt),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (upper) + mismatch')
subplot(6,1,6)
plot(t,th_noise_l(cnt,:) + th_mismatch_l(cnt),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (lower) + mismatch')

size(tex_rs_chunk)
size(trialNoise_filt_scaled)
size(sensorxtalk)
size(sensorxtalk_flicker)

tex_rs_chunk = tex_rs_chunk + trialNoise_filt_scaled + sensorxtalk + sensorxtalk_flicker;


[u,q,q_node,delta,vthp,vthn,y,LSB,spktimes_on,spktimes_off,spktimes_node] = LCSDQ(tex_rs_chunk,t,max_maxSR,N,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,th_noise_h,th_noise_l,th_mismatch_h,th_mismatch_l);
if wta == 1
    hem.spktimes = spktimes_node;
else
    hem.spktimes_on = spktimes_on;
    hem.spktimes_off = spktimes_off;
end

%% figures for paper - non-idealities
figure
subplot(4,1,1)
plot(t,sensorxtalk(1,:),'linewidth',2)
title('Signal conditioning output noise (white)')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(4,1,2)
plot(t,sensorxtalk_flicker(1,:),'linewidth',2)
title('Signal conditioning output noise (flicker)')
xlabel('Time (s)');
ylabel('Amplitude (V)');
subplot(4,1,3)
plot(t,LSB/2 + th_noise_h(1,:) + th_mismatch_h(1),'linewidth',2); hold on;
line([0 2], [LSB/2 LSB/2],'LineStyle','--','Color','r');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (upper) + mismatch')
xlim([0 2]);
legend('Actual Threshold', 'Ideal Threshold');
subplot(4,1,4)
plot(t,LSB/2 + th_noise_l(1,:) + th_mismatch_l(1),'linewidth',2); hold on;
line([0 2], [LSB/2 LSB/2],'LineStyle','--','Color','r');
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Threshold noise (lower) + mismatch')
xlim([0 2]);
legend('Actual Threshold', 'Ideal Threshold');

%% figures
ff2=figure(2);
subplot(3,1,1);
plot(t,y(1,:),'linewidth',2)
hold on
stairs(t,u(1,:),'linewidth',2);
stairs(t,vthp(1,:),'linewidth',2);
stairs(t,vthn(1,:),'linewidth',2);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('LC Quantization');
legend('Input signal','u','vthp','vthn');
%ylim([0 1]);
xlim([0 stim.tlength]);
subplot(3,1,2);
plot(t,delta(1,:),'linewidth',2)
xlabel('Time (s)');
ylabel('Amplitude (V)');
line([0 4], [LSB/2+th_mismatch_h(1) LSB/2+th_mismatch_h(1)],'LineStyle','--','Color','r');
line([0 4], [-LSB/2+th_mismatch_l(1) -LSB/2+th_mismatch_l(1)],'LineStyle','--','Color','r');
legend('delta','vth = LSB/2', 'vth = -LSB/2');
xlim([0 stim.tlength]);
ylim([-LSB LSB]);
subplot(3,1,3);
stem(t,q(1,:),'linewidth',2);
xlabel('Time (s)');
legend('q');
xlim([0 stim.tlength]);

if size(y,1) > 1
    ff12=figure(12);
    subplot(3,1,1);
    plot(t,y(2,:),'linewidth',2)
    hold on
    stairs(t,u(2,:),'linewidth',2);
    stairs(t,vthp(2,:),'linewidth',2);
    stairs(t,vthn(2,:),'linewidth',2);
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    title('LC Quantization');
    legend('Input signal','u','vthp','vthn');
    %ylim([0 1]);
    xlim([0 stim.tlength]);
    subplot(3,1,2);
    plot(t,delta(2,:),'linewidth',2)
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
    line([0 4], [LSB/2+th_mismatch_h(2) LSB/2+th_mismatch_h(2)],'LineStyle','--','Color','r');
    line([0 4], [-LSB/2+th_mismatch_l(2) -LSB/2+th_mismatch_l(2)],'LineStyle','--','Color','r');    
    legend('delta','vth = LSB/2', 'vth = -LSB/2');
    xlim([0 stim.tlength]);
    ylim([-LSB LSB]);
    subplot(3,1,3);
    stem(t,q(2,:),'linewidth',2);
    xlabel('Time (s)');
    legend('q');
    xlim([0 stim.tlength]);
end    

if crf==1
    cnt=hem_cnt;
else
    %if mr_cnt > 64
    %    cnt=64;  % fit as many heminodes
    %else
        cnt=mr_cnt;
    %end
end

%xdim=16; ydim=12;
xdim=sqrt(mr_cnt); ydim=sqrt(mr_cnt);
ff23=figure(23);
clf(ff23);
set(gcf,'position',[10,10,1800,1200])
for h=1:cnt
    %subplot(cnt+1,1,h)
    %subplot(sqrt(mr_cnt)+1,sqrt(mr_cnt),h);
    subplot(ydim+1,xdim,h);
    stem(t,q(h,:),'linewidth',2);  
    xlim([0 stim.tlength]);
    xlabel('Time (s)');
    title(['Sensor #',num2str(h)])
end
%subplot(cnt+1,1,cnt+1)
%subplot(sqrt(mr_cnt)+1,sqrt(mr_cnt),cnt+1)
%stem(t,q_node(:),'linewidth',2); 
%title('Spikes from node of Ranvier (WTA)')

ff24=figure(24);
clf(ff24);
set(gcf,'position',[10,10,1800,1200])
for h=1:cnt
    %subplot(cnt+1,1,h)
    %subplot(sqrt(mr_cnt),sqrt(mr_cnt),h);
    subplot(ydim,xdim,h);
    plot(t,y(h,:),'linewidth',2);  
    xlim([0 stim.tlength]);
    title(['Sensor #',num2str(h)])
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
end


for h=1:cnt
    SR{h} = gradient(y(h,:))/(t(2)-t(1));
    maxSR(h) = max(abs(SR{h}));   
    tlooprequired(h) = LSB/maxSR(h);
    disp(['Required tloop to avoid slope overload: <',num2str(tlooprequired(h))]);
end
%hem.tlooprequired = tlooprequired;
%hem.maxSR = maxSR;


% 
% ff3=figure(3);
% plot(t,SR{1},'linewidth',2)
% 
% ff4=figure(4);
% histogram(u);He

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
    spikecount_on(h) = numel(spktimes_on{h});
    spikecount_off(h) = numel(spktimes_off{h});
    disp(['for heminode # ',num2str(h),': ON spikerate = ',num2str(spikecount_on(h)/tlength),' Hz ',': OFF spikerate = ',num2str(spikecount_off(h)/tlength),' Hz']); 
end
hem.spikecount_on = spikecount_on;
hem.spikecount_off = spikecount_off;

%% display number of samples for all heminodes
for h=1:cnt
    %samplecnt(h) = length(find(q(h,:)==1));
    samplecnt_on(h) = numel(spktimes_on{h});
    samplecnt_off(h) = numel(spktimes_off{h});
    disp(['LCADC   N=',num2str(N),'  LSB= ',num2str(LSB),'V   texture:',tex_name{tsel},'  ON samples for hem #',num2str(h),' = ',num2str(samplecnt_on(h)),'  OFF samples for hem #',num2str(h),' = ',num2str(samplecnt_off(h))]);
end
hem.samplecnt_on = samplecnt_on;
hem.samplecnt_off = samplecnt_off;

if fileprint == 1
    % fname1 = (['figures/LCtexture-',tex_name{tsel},'-','s_fig1']);
     fname2 = ([current_loc,'figures/spikegen/LCtexture-',tex_name{tsel},'- N=',num2str(N),'s_fig2']);
    % fname3 = (['figures/LCtexture-',tex_name{tsel},'-','s_fig3']);
    % fname4 = (['figures/LCtexture-',tex_name{tsel},'-','s_fig4']);
    % print(ff1,'-dpng','-r600',fname1);
     print(ff2,'-dpng','-r600',fname2);
    % print(ff3,'-dpng','-r600',fname3);
    % print(ff4,'-dpng','-r600',fname4);
end
toc

end