function [maxSR,sig_ave,sig_pp] = print_textures(current_loc,stim,MR_type,plot_trials,hem_cnt)


tex_rs = stim.tex_rs;
tex_rsdiff = stim.tex_rsdiff;
tex_name = stim.tex_name;
L = stim.L;
t = stim.t;
ti = stim.ti;
tlength = stim.tlength;
fs_stim = stim.fs_stim;
stim_os = stim.stim_os;
%hem_cnt=16;


ff10 = figure(10);
set(gcf,'position',[10,10,1800,1200])

% print all texture for 1st heminode, 1st trial
for k = 1:size(tex_rs,1)

    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{k,1}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{k,1}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(3,4,k)
    subplot(ceil(sqrt(size(tex_rs,1))),ceil(sqrt(size(tex_rs,1))),k)
    plot(t,tex_rs{k}(1,:),'linewidth',2);
    title([tex_name{k}]);
    xlabel('Time (s)');
    %ylabel('Amplitude (V)');
    ylim([0 1]);
end

% %% PRINT ONLY HALF
% ff10 = figure(10);
% set(gcf,'position',[10,10,1800,1200])
% 
% % print all texture for 1st heminode
% for k = 1:size(tex_rs,1)/2
% 
%     if MR_type == 'FA'        
%         tex_rs_chunk = tex_rsdiff{k}(1,:); % print normalized sensor derivative
%     elseif MR_type == 'SA'
%         tex_rs_chunk = tex_rs{k}(1,:); % print normalized sensor amplitude  
%     else
%         tex_rs_chunk = 0; % get whole chunk
%     end  
%     
%     %subplot(3,4,k)
%     subplot(ceil(sqrt(size(tex_rs,1)/2)),ceil(sqrt(size(tex_rs,1)/2)),k)
%     plot(t,tex_rs{k}(1,:),'linewidth',2);
%     title([tex_name{k}]);
%     xlabel('Time (s)');
%     %ylabel('Amplitude (V)');
% end

ff19 = figure(19);
set(gcf,'position',[10,10,1800,1200])
% print all heminodes for 1st texture, 1st trial
for h = 1:hem_cnt

    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{k,1}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{k,1}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(4,4,h)
    subplot(ceil(sqrt(hem_cnt)),ceil(sqrt(hem_cnt)),h)
    plot(t,tex_rs{1}(h,:),'linewidth',2);
    title([tex_name{1}]);
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
end 

if plot_trials==1
    trialrange=size(tex_rs,2);
else
    trialrange=40;
end

trialrange
ff20 = figure(20);
set(gcf,'position',[10,10,1800,1200])
% print all trials from 1st texture (1st heminode)

for trial = 1:trialrange

    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{1,trial}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{1,trial}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(4,4,h)
    subplot(ceil(sqrt(size(tex_rs,2))),ceil(sqrt(size(tex_rs,2))),trial)
    plot(t,tex_rs{1,trial}(1,:),'linewidth',2);
    title([tex_name{1},' trial # ',num2str(trial)]);
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
end 

ff21 = figure(21);
set(gcf,'position',[10,10,1800,1200])
% print all trials from 1st texture (1st heminode)
for trial = 1:trialrange

    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{size(tex_rs,1),trial}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{size(tex_rs,1),trial}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(4,4,h)
    subplot(ceil(sqrt(size(tex_rs,2))),ceil(sqrt(size(tex_rs,2))),trial)
    plot(t,tex_rs{size(tex_rs,1),trial}(1,:),'linewidth',2);
    title([tex_name{size(tex_rs,1)},' trial # ',num2str(trial)]);
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
end 

ff22 = figure(22);
set(gcf,'position',[10,10,1800,1200])
% print all trials from 1st texture (1st heminode)
for trial = 1:trialrange

    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{size(tex_rs,1)-1,trial}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{size(tex_rs,1)-1,trial}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(4,4,h)
    subplot(ceil(sqrt(size(tex_rs,2))),ceil(sqrt(size(tex_rs,2))),trial)
    plot(t,tex_rs{size(tex_rs,1)-2,trial}(1,:),'linewidth',2);
    title([tex_name{size(tex_rs,1)-2},' trial # ',num2str(trial)]);
    xlabel('Time (s)');
    ylabel('Amplitude (V)');
end 


%% plot spectrum of all textures, 1st heminode
ff11=figure(11);
set(gcf,'position',[10,10,1800,1200])

for k = 1:size(tex_rs,1)
    
    if MR_type == 'FA'        
        tex_rs_chunk = tex_rsdiff{k,1}(1,:); % print normalized sensor derivative
    elseif MR_type == 'SA'
        tex_rs_chunk = tex_rs{k,1}(1,:); % print normalized sensor amplitude  
    else
        tex_rs_chunk = 0; % get whole chunk
    end  
    
    %subplot(3,4,k);
    subplot(ceil(sqrt(size(tex_rs,1))),ceil(sqrt(size(tex_rs,1))),k)
    Fv = linspace(0,1,fix(length(tex_rs_chunk)/2)+1)*(fs_stim*stim_os/2);
    Iv = 1:length(Fv);
    FFTrs = fft(tex_rs_chunk)/length(tex_rs_chunk);
    plot(10*log10(abs(FFTrs(Iv))*2));
    xlabel('Frequency (Hz)');
    xlim([0 10e3]);
    ylabel('Magnitude (dB)');
    title(['Spectrum of texture: ',tex_name{k}]);
end

ff12=figure(12);   % plot 1 texture spectrum only
set(gcf,'position',[10,10,1800,1200])


k=1;  %% select texture to plot
if MR_type == 'FA'        
    tex_rs_chunk = tex_rsdiff{k,1}(1,:); % print normalized sensor derivative
elseif MR_type == 'SA'
    tex_rs_chunk = tex_rs{k,1}(1,:); % print normalized sensor amplitude  
else
    tex_rs_chunk = 0; % get whole chunk
end  

%subplot(3,4,k);
n=length(tex_rs_chunk);
%FFTrs = fft(tex_rs_chunk)/length(tex_rs_chunk);
FFTrs = fft(tex_rs_chunk);
PSD=FFTrs.*conj(FFTrs)/n;
freq=n/(n)*(0:n);
L=1:floor(n/2);
plot(freq(L),PSD(L));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(['Spectrum of texture: ',tex_name{k}]);


% ff14=figure(14);   % plot 1 texture spectrum only
% set(gcf,'position',[10,10,1800,1200])
% 
% 
% k=50;  %% select texture to plot
% if MR_type == 'FA'
%     tex_rs_chunk = tex_rsdiff{k}(1,:); % get whole chunk
% elseif MR_type == 'SA'
%     tex_rs_chunk = tex_rs{k}(1,:); % get whole chunk   
% else
%     tex_rs_chunk = 0; % get whole chunk
% end 

%subplot(3,4,k);
n=length(tex_rs_chunk);
%FFTrs = fft(tex_rs_chunk)/length(tex_rs_chunk);
FFTrs = fft(tex_rs_chunk);
PSD=FFTrs.*conj(FFTrs)/n;
freq=n/(n)*(0:n);
L=1:floor(n/2);
plot(freq(L),PSD(L));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(['Spectrum of texture: ',tex_name{k}]);

ff13=figure(13);
set(gcf,'position',[10,10,1800,1200])

for k = 1:size(tex_rs,1)
    for trial=1:size(tex_rs,2)
        for h=1:size(tex_rs{k,trial},1)  % get number of rows == number of heminodes
            if MR_type == 'FA'
                tex_rs_chunk = tex_rsdiff{k,trial}(h,:); % get whole chunk
            elseif MR_type == 'SA'
                tex_rs_chunk = tex_rs{k,trial}(h,:); % get whole chunk   
            else
                tex_rs_chunk = 0; % get whole chunk
            end   

            SR{h,k,trial} = gradient(tex_rs_chunk)/(t(2)-t(1));
            maxSR(h,k,trial) = max(abs(SR{h,k,trial}));
            sig_ave(h,k,trial) = mean(tex_rs_chunk);
            start_i=round(length(tex_rs_chunk)*0.1) % start at 0.2s to discard ramp up from 0
            sig_pp(h,k,trial) = max(tex_rs_chunk(start_i:end)) - min(tex_rs_chunk(start_i:end));
        end
    end
    %subplot(3,4,k);
    subplot(ceil(sqrt(size(tex_rs,1))),ceil(sqrt(size(tex_rs,1))),k)
    plot(t,SR{1,k},'linewidth',2)
    title(['dV/dt of 1st heminode: ',tex_name{k}]);
end



end