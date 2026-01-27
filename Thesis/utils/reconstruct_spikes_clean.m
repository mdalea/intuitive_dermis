%% Author: Mark Daniel Alea
% mark.alea@kuleuven.be
% 04/02/2023
% reads .bs2 files and display as spatiotemporal spikes


folderPath = pwd; % replace with path to dataset
name = 'N_3b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_180mvpp_250mvofs';
%name = 'N_4b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_180mvpp_250mvofs_onedac';
%name = 'N_3b_allchannels_ALLCLASS_auto_nospkfilter_cleaned_12x16_vibration'
sample=1;  % choose which sample to plot
chunks=0;  % 0 if dataset without the -x
%name = 'latiti'
%name = 'N_3b_allchannels_1p80mvpp_250mvofs_cl1_tr1';

if chunks==0
    append = '';
else
    append = ['-',num2str(chunks)];
end

filename = [folderPath,'/out/',name,'/',num2str(sample),append,'.bs2']

TD = Read_Ndataset(filename);

tduration=0.4; % 0.4s for texture; 2s for vibration
t=[0:1:tduration*1e6]; % 400ms 1us period
sig = zeros(192,length(t)); % 1us sampling period to reflect spiketime resolution
N=3;
delta_hi=1/2^N;
delta_lo=2/2^N;

figure
hold all;
for spikeCount = 1:length(TD.ts)
    if TD.p(spikeCount)==1  % ON spike
        %xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)
        %sig(xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+1,TD.ts) = sig(xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+1,TD.ts) + delta;
        plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
    else   % OFF spike
       %sig(xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+1,TD.ts) = sig(xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+1,TD.ts) - delta;       
       plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
            [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
    end 
end

grid on
xlabel('Time (ms)');
ylabel('Taxel Number');
ylim([0 191])
title('INTUITIVE spike output')

micasplot_large

for x=1:16
    for y=1:12
        TD_t = [];
        TD_t.x = TD.x(TD.x==x-1 & TD.y==y-1)
        TD_t.y = TD.y(TD.x==x-1 & TD.y==y-1)
        TD_t.p = TD.p(TD.x==x-1 & TD.y==y-1)
        TD_t.ts = TD.ts(TD.x==x-1 & TD.y==y-1)
        for spikeCount = 1:length(TD_t.ts)
            length(TD_t.ts)
            
            if TD_t.p(spikeCount)==1  % ON spike
                %xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)
                sig(xytolin_customxy(TD_t.x(spikeCount),TD_t.y(spikeCount),16)+1,TD_t.ts(spikeCount):length(t)) = sig(xytolin_customxy(TD_t.x(spikeCount),TD_t.y(spikeCount),16)+1,TD_t.ts(spikeCount):length(t)) + delta_hi;
                %plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
                %    [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k');
            else   % OFF spike
               sig(xytolin_customxy(TD_t.x(spikeCount),TD_t.y(spikeCount),16)+1,TD_t.ts(spikeCount):length(t)) = sig(xytolin_customxy(TD_t.x(spikeCount),TD_t.y(spikeCount),16)+1,TD_t.ts(spikeCount):length(t)) - delta_lo;       
               %plot([TD.ts(spikeCount)/1e3 TD.ts(spikeCount)/1e3], ...
               %     [xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)-0.4 xytolin_customxy(TD.x(spikeCount),TD.y(spikeCount),16)+0.4], 'k','Color','red');
            end 
        end
    end
end

figure
%hold on;
for i=1:192
    subplot(12,16,i)
    plot(t,sig(i,:))
    maxsig(i)=max(sig(i,:))
    minsig(i)=min(sig(i,:))
end

    voutpp=maxsig-minsig
    max(voutpp)

%micasplot



