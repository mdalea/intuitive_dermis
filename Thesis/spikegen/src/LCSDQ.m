%% function

function [u,q,q_node,delta,vthp,vthn,y_in,LSB,spktimes_on,spktimes_off,spktimes_node] = LCSDQ(y,t,max_maxSR,N,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,th_noise_h,th_noise_l,th_mismatch_h,th_mismatch_l)
%% what does this code do?
%as per basic equations,models a level-crossing multibit DAC ADC

if crf == 1
    LSB = max(mr_in_hem)*(1/(2^N));  % higher threshold to normalize firing rate    
    tloop = LSB / max_maxSR
else
    LSB = 1/(2^N);
    tloop = LSB / max_maxSR
end

timestep = t(2) - t(1);  % assume uniformly spaced time
% set time granularity as tloop
if tloop < timestep
    tloopnew = timestep;
    disp(['tloop = ',num2str(tloop),' < timestep = ',num2str(timestep),' setting tloop = timestep']);
else
    if round(tloop/timestep) ~= tloop/timestep
        tloopnew = round(tloop/timestep)*timestep;
        disp(['timestep is ',num2str(timestep),' Rounding off tloop to ',num2str(tloopnew)]);
    else
        tloopnew = tloop;
        disp(['Keeping tloop to ',num2str(tloopnew)]); 
    end
end

kmod = round(tloopnew/timestep);

if crf==1
    y_sum = zeros(hem_cnt,length(t));
    for h=1:hem_cnt
        for i=1:numel(mr_index{h})
           y_sum(h,:) = y_sum(h,:) + y(mr_index{h}(i),:);
        end
        figure(6)
        %plot(DT*(1:round(gTmax/DT)),g_syn);
        h(h) = plot(t,y_sum(h,:)); hold on; 
        ylabel('Input to heminode (fullscale units)');
        xlabel('Time (s)');
    end
    %lgd = legend(h,'Hem #1','Hem #2','Hem #3','Hem #4');
    %lgd.FontSize = 9; lgd.Location = 'southwest';
    hold off;
end

if crf==1
    cnt=hem_cnt;
else
    cnt=mr_cnt;
end

size(th_noise_h)
size(th_mismatch_h)


%% code logic
q = zeros(cnt,length(t));
q_node = zeros(1,length(t));
u = zeros(cnt,length(t)); %feedback DAC value
vthp = zeros(cnt,length(t));
vthn = zeros(cnt,length(t));
delta = zeros(cnt,length(t));
%dir = zeros(cnt,length(t)); % direction value
vth0 = LSB/2;
u(:,1) = vth0; % set DAC to be at the middle of vthp and vthn by 1 LSB
%u(:,1) = 0;
%u(:,1) = y(:,1);   % prevent slope overload at the start
%vthp(:,1) = u(:,1) + vth0;
%vthn(:,1) = u(:,1) - vth0;

y_in = zeros(cnt,length(t));
for h=1:cnt
    if crf==1
        y_in(h,:) = y_sum(h,:);
    else
        y_in(h,:) = y(h,:);
    end
end


%recursive equations for SDQ
%% no WTA implementation yet
for k = 2:length(t)
    for h=1:cnt
        delta(h,k) = y_in(h,k-1) - u(h,k-1);
        u(h,k) = u(h,k-1);
        %vthp(h,k) = u(h,k-1) + vth0;
        %vthn(h,k) = u(h,k-1) - vth0;
        q(h,k) = 0;
        dir(h,k) = 0;
        if mod(k,kmod) == 0
            if delta(h,k) >= vth0 + th_noise_h(h,k) + th_mismatch_h(h)
                q(h,k) = 1;
                %dir(h,k) = 1;
                u(h,k) = u(h,k) + vth0;
            elseif delta(h,k) < -vth0 + th_noise_l(h,k) + th_mismatch_l(h)
                q(h,k) = -1;
                %dir(h,k) = -1;     
                u(h,k) = u(h,k) - vth0;
            end    
        end
    end
    
end  

spktimes_on = [];
spktimes_off = [];
%spktimes_dir = [];
for h=1:cnt   
    % get times where there is a spike, denormalize to seconds
    %DT = 1e3/0.0001;
    DT = 1e3/0.000001;
    %spktimes_dir{h} = dir(h,find(q(h,:)==1));

    spktimes_on{h} = t(find(q(h,:)==1));
    spktimes_on{h} = single(spktimes_on{h});
    spktimes_on{h} = round(spktimes_on{h}*DT)/DT; % round up to 4 decimal places for classifier
    %spktimes_dir{h} = spktimes_dir{h}(spktimes{h}>0.1);  % discard spikes during DAC settling
    spktimes_on{h} = spktimes_on{h}(spktimes_on{h}>0.01);  % discard spikes during DAC settling

    spktimes_off{h} = t(find(q(h,:)==-1));
    spktimes_off{h} = single(spktimes_off{h});
    spktimes_off{h} = round(spktimes_off{h}*DT)/DT; % round up to 4 decimal places for classifier
    %spktimes_dir{h} = spktimes_dir{h}(spktimes{h}>0.1);  % discard spikes during DAC settling
    spktimes_off{h} = spktimes_off{h}(spktimes_off{h}>0.01);  % discard spikes during DAC settling
end

spktimes_node = [];
spktimes_node_dir = [];

%% no WTA implementation yet
%spktimes_node_dir{1} = dir(1,find(q_node(:)==1));

% spktimes_node{1} = t(find(q_node(:)==1));
% spktimes_node{1} = single(spktimes_node{1});
% spktimes_node{1} = round(spktimes_node{1}*DT)/DT; % round up to 4 decimal places for classifier
% spktimes_node_dir{1} = spktimes_node_dir{1}(spktimes_node{1}>0.1);  % discard spikes during DAC settling
% spktimes_node{1} = spktimes_node{1}(spktimes_node{1}>0.1);


end 