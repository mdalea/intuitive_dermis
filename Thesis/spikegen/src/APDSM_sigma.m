%% function

function [u,q,q_node,delta,vth,g_syn,y_in,spktimes,spktimes_node] = APDSM_sigma(y,t,max_maxSR,N,DT,LSB,LSB_th,synfilt_en,tau_syn,tau_refr,mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,th_noise,th_mismatch,tau_refr_mismatch)
%% what does this code do?
%as per basic equations,models an APDSM in [Yoon,IEEETNNL17]

if crf == 1
    LSB_nom = max(mr_in_hem)*(1/(2^N));  % higher threshold to normalize firing rate    
    trefract = (LSB_nom*tau_refr) / (LSB_nom + tau_refr*max_maxSR) / 3
else
    LSB_nom = 1/(2^N);
    trefract = (LSB_nom*tau_refr) / (LSB_nom + tau_refr*max_maxSR) / 3
end

timestep = t(2) - t(1);  % assume uniformly spaced time
% set time granularity as tloop
if trefract < timestep
    trefractnew = timestep;
    disp(['trefract = ',num2str(trefract),' < timestep = ',num2str(timestep),' setting trefract = timestep']);
else
    if round(trefract/timestep) ~= trefract/timestep
        trefractnew = round(trefract/timestep)*timestep;
        disp(['timestep is ',num2str(timestep),' Rounding off trefract to ',num2str(trefractnew)]);
    else
        trefractnew = trefract;
        disp(['Keeping trefract to ',num2str(trefractnew)]);        
    end
end

kmod = round(trefractnew/timestep);

tau_syn
%% exponential kernel for synaptic kernel
gTmax = 5*tau_syn;
%gVar_inp_min = LSB/2;

% Creating exponential kernel
TempT = zeros(1,length(t)); % Initialization
g_syn = zeros(1,length(t)); % Initialization

% Initial Parameters
g_syn(1) = LSB_nom/2;
TempT(1)=0;

% Template
for i=2:length(t)
    TempT(i)=TempT(i-1)+DT;
    %if abs(TempT(i)-1)<DT/100; 
    %    g_syn(i-1)=gVar_inp_min; 
    %end
    g_syn(i)= g_syn(i-1) - DT/tau_syn * g_syn(i-1);
end


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

%% exponential kernel for refractory kernel
for h=1:cnt
    gTmax(h) = 5*(tau_refr + tau_refr_mismatch(h));
end
%gVar_inp_min = LSB/2;

% Creating exponential kernel
TempT = zeros(h,length(t)); % Initialization
g_refr = zeros(h,length(t)); % Initialization

% Initial Parameters
for h=1:cnt
    g_refr(h,1) = LSB(h)/2;
    TempT(h,1)=0;
end

% Template
for h=1:cnt
    for i=2:length(t)
        TempT(h,i)=TempT(h,i-1)+DT;
        %if abs(TempT(i)-1)<DT/100; 
        %    g_syn(i-1)=gVar_inp_min; 
        %end
        g_refr(h,i)= g_refr(h,i-1) - DT/(tau_refr+tau_refr_mismatch(h)) * g_refr(h,i-1);
    end
end


figure(5)
%plot(DT*(1:round(gTmax/DT)),g_syn);
plot(t,g_syn); hold on;
plot(t,g_refr(1,:));
plot(t,g_refr(2,:));
plot(t,y(1,:));
%DS_lim = round(gTmax/DT);
DS_lim = length(t);

%% code logic
q = zeros(cnt,length(t)+DS_lim);
q_node = zeros(1,length(t)+DS_lim);
u = zeros(cnt,length(t)+DS_lim); %feedback DAC value
vth = zeros(cnt,length(t)+DS_lim);
delta = zeros(cnt,length(t)+DS_lim);
u(:,1) = 0; %speed up start-up

% init_steps(:) = round(y(:,1)/max(g_syn))   % compute how many steps of kernel at the beginning (fast start-up, no slope overload)
% u(:,1:1+DS_lim-1) = init_steps(:)*g_syn(1:end);  % prevent slope overload at the start

% %% TEST STEP RESPONSE (COMMENT OUT)
% step = zeros(1,length(t));
% step(length(step)/10:end) = 1;
% y(:,:) = zeros(cnt,length(t));
% for h=1:cnt
%     y(h,:) = step;
% end
% %%


y_in = zeros(cnt,length(t));
for h=1:cnt
    if crf==1
        y_in(h,:) = y_sum(h,:);
    else
        y_in(h,:) = y(h,:);
    end
end

if (synfilt_en==1)
    %% CONVOLVE INPUT WITH SYNAPTIC KERNEL (LPF)
    for h=1:cnt
        if round(gTmax/DT) > length(y_in(h,:))   % tau too large, tail exceeds sensor recording time
            y_sigma = conv(y_in(h,:),g_syn(1:length(y_in(h,:)))); 
        else
            y_sigma = conv(y_in(h,:),g_syn(1:round(gTmax/DT))); 
        end
        scale_factor = max(y_in(h,:)) / max(y_sigma(1:length(y_in(h,:))));
        y_in(h,:) = scale_factor.*y_sigma(1:length(y_in(h,:)));
    end
end

% assumes total loop delay is 1/(os_stim*fs_stim) -> comparison is done
% every sample
for k=2:length(t)
    crossing=0;  % indicates if one heminode spikes
    for h=1:cnt
        delta(h,k) = y_in(h,k) - u(h,k);
        vth(h,k) = u(h,k) + LSB_nom/2;
        %if y(k) > vth(k-1)
        if mod(k,kmod) == 0
            %if y(h,k) > vth(h,k)
            if delta(h,k) > LSB_th/2 + th_noise(h,k) + th_mismatch(h) % use nom LSB so mismatch does not pile
                q(h,k) = 1;
                if wta == 0
                    u(h,k:k+DS_lim-1) = u(h,k:k+DS_lim-1) + g_refr(h,1:end);  % add refractory kernels, not synaptic 
                else
                    crossing=1;
                end
            end    
        end
    end
    
    if (wta == 1) && (crossing==1)
        q_node(k) = 1;
        u(:,k:k+DS_lim-1) = u(:,k:k+DS_lim-1) + g_refr(h,1:end);
    end
   
end 

spktimes = [];
for h=1:cnt   
    % get times where there is a spike, denormalize to seconds
    %DT = 1e3/0.0001;
    DT = 1e3/0.000001;
    spktimes{h} = t(find(q(h,:)==1));
    spktimes{h} = single(spktimes{h});
    spktimes{h} = round(spktimes{h}*DT)/DT; % round up to 4 decimal places for classifier
    spktimes{h} = spktimes{h}(spktimes{h}>0.1);     % discard spikes during DAC settling  
end

spktimes_node = [];
spktimes_node{1} = t(find(q_node(:)==1));
spktimes_node{1} = single(spktimes_node{1});
spktimes_node{1} = round(spktimes_node{1}*DT)/DT; % round up to 4 decimal places for classifier
spktimes_node{1} = spktimes_node{1}(spktimes_node{1}>0.1);

% truncate signals
delta = delta(:,1:length(t));
vth = vth(:,1:length(t));
u = u(:,1:length(t));
q = q(:,1:length(t));
q_node = q_node(1:length(t));

end