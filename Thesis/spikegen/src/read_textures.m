%% handles exponential TF only for LMT & LMT (filt)
%% handles reduced classes only for LMT(Filt)
function [stim] = read_textures(file_loc,dataset,tex_cnt,trialcnt,randtex,randtex_indices,fs_stim,stim_os,mr_cnt,scanv,scandir,exp_tf,spatialres)

    K=10;  % no. of top classes to use
    %K_ind = [29 10 30 38 39 40 26 42 45 17 16]; % indices, get from topK results. Change as you please
    %K_ind = [56 68 45 44 63 52 38 16 40 22]; %[56 45 68 17  7 40 20 19 52 31];
    K_ind = [14  5 55 58 15 44 52 40  8 10];
    K_ind = K_ind + 1;  % index 0 start
    %spatialres=0.2; %mm

if dataset == "Kursun"
    
    %% from Kursun tactile dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % assume csv file is read already into X matrix
    fs_stim = 200; % original signal frequency
    stim_os = 25;  % match overall LMT sampling rate (5kHz)
    
    X = csvread(file_loc);

    %DT=1/(fs_stim*stim_os);
    tlength_orig=20;   % original time length
    %tdelay=1;  % phase shift between seen sensor values by heminode
    scanv_orig = 50; % scanning velocity from Kursun paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;   % adjust duration of recording dependent on scanning speed (i.e., longer if slower, shorter if faster)
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    %tdelaycnt = 16; % no. of delayed samples per texture & trial
    tdelaycnt = 2*sqrt(mr_cnt);
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 12;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 1;   
    
        
    for k = 1:size(X,1)
       disp(['Reading stimulus # ',num2str(k)]);
       tex{k,1} = X(k,:); % csv file divides different textures into different row
       
       L = length(tex{1});
       t = (0:L-1)*(1/fs_stim);  % get duration of stimulus recording - expect 20seconds, 5 drum rotations
       %[tex{k,1},~] = lsim(pos_filt,tex{k,1},t);  
    end
    

    tex_name = {'fabric-1','aluminum film','fabric-2','fabric-3','moquette-1','moquette-2','fabric-4','sticky fabric-5','sticky fabric','sparkle paper-1','sparkle paper-2','toy-tire rubber'};    
    
    for k = 1:numel(tex_name)
        rs = resample(tex{k,1},stim_os*length(tex{k,1}),length(tex{k,1}));   % oversample recording
        tex{k,1} = rs;   
        texmax(k) = max(rs);
        texmin(k) = min(rs);            
    end
    
    texmaxmax = max(max(texmax))
    texminmin = min(min(texmin))    
    
    for k = 1:numel(tex_name)
        rs = tex{k,1};
        
        if exp_tf==0
            rs = rs + abs(texminmin);   % push negative values above zero
            rs = (1/(texmaxmax-texminmin)).*rs;                
        else
            rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
            rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
        end
           
        tex{k,1} = rs;
        texmax(k) = max(rs);
        texmin(k) = min(rs);    
    end
    
    texmaxmax = max(max(texmax))
    texminmin = min(min(texmin))     
     
    for k = 1:numel(tex_name)     
       rs = tex{k,1};
       if exp_tf==1   % normalize after exp()
            rs = (1/(texmaxmax-texminmin)).*rs;
       end    
        
       %tex_rs{k,trial} = rs;
       temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
       temp = temp';
       %tlength_des=18;  % new desired time length
       tlength_des = 0.833*tlength;
       for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
            %disp(['--->Creating time chunk # ',num2str(t)]);  
            temp_t = [];
            for i=1:round(tlength_des/tdelay)
                temp_t = [temp_t,temp(i+t-1,:)];
            end
            tex_rs_t{k,1,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
       end
    end

    size(tex_rs_t)

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir); 
     
    %tsel = 1; % select which texture
    %n = 1; % divide recording into 1second bites; which 1s bite to use
    %L = length(tex_rs{tsel}(n,:));
    L = length(tex_rs_t{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector
    
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    
    
    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling

elseif dataset == "Kursun_filt"
    
    %% from Kursun tactile dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % assume csv file is read already into X matrix
    fs_stim = 200; % original signal frequency
    stim_os = 25;  % match overall LMT sampling rate (5kHz)
    
    X = csvread(file_loc);

    %DT=1/(fs_stim*stim_os);
    tlength_orig=20;   % original time length
    %tdelay=1;  % phase shift between seen sensor values by heminode
    scanv_orig = 50; % scanning velocity from Kursun paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;   % adjust duration of recording dependent on scanning speed (i.e., longer if slower, shorter if faster)
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    %tdelaycnt = 16; % no. of delayed samples per texture & trial
    tdelaycnt = 2*sqrt(mr_cnt);
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 12;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 1;   
    
    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
        
    for k = 1:size(X,1)
       disp(['Reading stimulus # ',num2str(k)]);
       tex{k,1} = X(k,:); % csv file divides different textures into different row
       
       L = length(tex{1});
       t = (0:L-1)*(1/fs_stim);  % get duration of stimulus recording - expect 20seconds, 5 drum rotations
       [tex{k,1},~] = lsim(pos_filt,tex{k,1},t);  
    end
    

    tex_name = {'fabric-1','aluminum film','fabric-2','fabric-3','moquette-1','moquette-2','fabric-4','sticky fabric-5','sticky fabric','sparkle paper-1','sparkle paper-2','toy-tire rubber'};    
    
    for k = 1:numel(tex_name)
        rs = resample(tex{k,1},stim_os*length(tex{k,1}),length(tex{k,1}));   % oversample recording
        tex{k,1} = rs;   
        texmax(k) = max(rs);
        texmin(k) = min(rs);            
    end
    
    texmaxmax = max(max(texmax))
    texminmin = min(min(texmin))    
    
    for k = 1:numel(tex_name)
        rs = tex{k,1};
        
        if exp_tf==0
            rs = rs + abs(texminmin);   % push negative values above zero
            rs = (1/(texmaxmax-texminmin)).*rs;                
        else
            rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
            rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
        end
           
        tex{k,1} = rs;
        texmax(k) = max(rs);
        texmin(k) = min(rs);    
    end
    
    texmaxmax = max(max(texmax))
    texminmin = min(min(texmin))     
     
    for k = 1:numel(tex_name)     
       rs = tex{k,1};
       if exp_tf==1   % normalize after exp()
            rs = (1/(texmaxmax-texminmin)).*rs;
       end    
        
       %tex_rs{k,trial} = rs;
       temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
       temp = temp';
       %tlength_des=18;  % new desired time length
       tlength_des = 0.833*tlength;
       for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
            %disp(['--->Creating time chunk # ',num2str(t)]);  
            temp_t = [];
            for i=1:round(tlength_des/tdelay)
                temp_t = [temp_t,temp(i+t-1,:)];
            end
            tex_rs_t{k,1,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
       end
    end

    size(tex_rs_t)

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir); 
     
    %tsel = 1; % select which texture
    %n = 1; % divide recording into 1second bites; which 1s bite to use
    %L = length(tex_rs{tsel}(n,:));
    L = length(tex_rs_t{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector
    
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    
    
    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling
elseif dataset == "Kursun_filt_timediv"
    %% divide 20second recording to 10 2s trials, just like what is used in Kursun paper
    
    %% from Kursun tactile dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % assume csv file is read already into X matrix
    fs_stim = 200; % original signal frequency
    stim_os = 25;  % match overall LMT sampling rate (5kHz)
    
    X = csvread(file_loc);

    %DT=1/(fs_stim*stim_os);
    tlength_orig=20;   % original time length (per trial)
    %tdelay=1;  % phase shift between seen sensor values by heminode
    scanv_orig = 50; % scanning velocity from Kursun paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;   % adjust duration of recording dependent on scanning speed (i.e., longer if slower, shorter if faster)
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10;  % divide into 10 chunks
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    %tdelaycnt = 16; % no. of delayed samples per texture & trial
    tdelaycnt = 2*sqrt(mr_cnt);
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 12;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 1;   

    tex_name = {'fabric-1','aluminum film','fabric-2','fabric-3','moquette-1','moquette-2','fabric-4','sticky fabric-5','sticky fabric','sparkle paper-1','sparkle paper-2','toy-tire rubber'};  
    
    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    trial_len = size(X,2)/10
    for k = 1:size(X,1)
       for trial = 1:10   % get 10 trials from 20s reconding (2s each)
        disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},' trial ',num2str(trial)]);
        tex{k,trial} = X(k,trial_len*(trial-1)+1:trial_len*(trial)); % csv file divides different textures into different row
 
         %% FILTER
        L = length(tex{1,1});  % assume all recording lengths are uniform
        t = (0:L-1)*(1/(fs_stim)); % time vector % get duration of stimulus recording -  20seconds is 5 drum rotations
        [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);        
       end
    end
    

  
    
    

    %%-----
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
           if exp_tf==1
              rs = exp(100*rs); 
           end  
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0    % no more negative values after the exp()
                if(texminmin<0)
                    rs = rs + abs(texminmin);   % push negative values above zero
                end
           else
                rs = rs - abs(texminmin);   % push lowest value to 0
           end

%            max(rs)
%            min(rs)
           %rs = 0.2 + 0.25.*(1/(texmaxmax-texminmin)).*rs;  % scale input signal to POSFET swing (0.25V), add POSFET offset (0.2V) due to mirror  
           rs = (1/(texmaxmax-texminmin)).*rs;

           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end

    %%----

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir); 
     
    %tsel = 1; % select which texture
    %n = 1; % divide recording into 1second bites; which 1s bite to use
    %L = length(tex_rs{tsel}(n,:));
    L = length(tex_rs_t{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector
    
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    
    
    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling

elseif dataset == "edge"
    
    fs_stim=5000;  % sampling rate in Saal model
    stim_os=240;   % to match oversampling rate in Kursun dataset

    %DT=1/(fs_stim*stim_os);
    tlength=1;   % original time length

    tex{1} = readmatrix([file_loc,'stim_stat_0.txt']);
    tex{2} = readmatrix([file_loc,'stim_stat_45.txt']);
    tex{3} = readmatrix([file_loc,'stim_stat_90.txt']);
    tex{4} = readmatrix([file_loc,'stim_stat_135.txt']);
    tex{5} = readmatrix([file_loc,'stim_stat_180.txt']);
    tex{6} = readmatrix([file_loc,'stim_stat_225.txt']);
    tex{7} = readmatrix([file_loc,'stim_stat_270.txt']);
    tex{8} = readmatrix([file_loc,'stim_stat_315.txt']);    
    size(tex{1})
    
    %% normalize 
    for k=1:numel(tex)
        for h=1:mr_cnt
            max(tex{k}(h,:));
            texmax(h,k) = max(tex{k}(h,:));
        end
    end
    texmaxmax = max(max(texmax))
    for k=1:numel(tex)
        tex{k} = tex{k}./texmaxmax;
    end
    
    for k=1:length(tex)
        for h=1:mr_cnt
            % if pressure is greater than some threshold, round it off
            % to max value -> ignore skin mechanics model
            if max(tex{k}(h,:)) > 0.1
                disp(['Resampling stimulus # ',num2str(k),' hem ',num2str(h)]);
                rs{k}(h,:) = resample((tex{k}(h,:))',stim_os*length(tex{k}(h,:)),length(tex{k}(h,:)));
                rs{k}(h,:) = (rs{k}(h,:))';
                rs{k}(h,:) = 0.75.*rs{k}(h,:)./max(tex{k}(h,:));   % normalize to 0.5 max
            else
                disp(['Ignoring stimulus # ',num2str(k),' hem ',num2str(h)]);
                rs{k}(h,:) = zeros(1,stim_os*length(tex{k}(h,:)));
            end
            size(rs{k})
        end
    end




    %% load stimulus data to all receptors
    for k=1:length(rs)
        for h=1:mr_cnt
            disp(['--> Filling stimulus # ',num2str(k),' hem ',num2str(h)]);
            tex_rs{k}(h,:) = rs{k}(h,:);                    
        end
    end



    %tsel = 1; % select which texture
    %n = 1; % divide recording into 1second bites; which 1s bite to use
    %L = length(tex_rs{tsel}(n,:));
    L = length(tex_rs{1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector

    % differentiate & rectify sensor voltage
    for k = 1:length(tex_rs)
        disp(['Reading stimulus # ',num2str(k)]);
       for h=1:mr_cnt
            disp(['--->Filling heminode # ',num2str(h),' (diff)']);
            temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
            pinds_diff = find(temp > 0);
            ninds_diff = find(temp < 0);  

            tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
            tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   

            tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
            tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
       end
    end

    tex_name = {'0 deg','45 deg','90 deg','135 deg','180 deg','225 deg','270 deg','315 deg'};

    
    stim.tex_rs = tex_rs;
    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    
    
elseif dataset == "LMT"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 108;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 1;
    
    tex_name = {'G1EpoxyRasterPlate','G1FineAluminumMesh','G1IsolatingFoilMesh','G1MeshFloorCloth', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2','G1ProfiledAluminumPlate','G1RhombAluminumMesh', ...
              'G1RhombAluminumMeshVersion2','G1RoundAluminumMesh','G1RubberMesh','G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta', 'G3AcrylicGlass', 'G3AluminumPlate', ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3EpoxyPlate', 'G3Glass', 'G3StainlessSteel', 'G4Bamboo', ...
              'G4Beech_0', 'G4CherryTree', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4CoreBeech', 'G4Cork', 'G4LaminatedWood', 'G4Larch', ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5RubberPlateVersion2', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2', 'G6IsolatingFoilVersion3', 'G6RedSleaze', ...
              'G6ScouringPadVersion1', 'G6ScouringPadVersion2', 'G6SheepSkin', 'G6Steelwool', ...
              'G6Towel', 'G6VelcroHooks', 'G6VelcroLoops', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7EpdmFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', 'G7MediumFoamVersion2', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8CarbonFoil', 'G8Cardboard', 'G8CardboardVersion2', 'G8GlitterPaperVersion1', ...
              'G8GlitterPaperVersion2', 'G8HighDensityPolyethylen', 'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2', 'G8PlayingCardPaper', 'G8RoughPaper', 'G8Tarp', ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9BlueSatin', 'G9FloorCloth', ...
              'G9FloorClothVersion2', 'G9FloorClothVersion3', 'G9GreenVelvet', 'G9Jeans', ...
              'G9Kashmir', 'G9Leather', 'G9LeatherVersion2', 'G9RedFleeze', ...
              'G9RedVelvet', 'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3', 'G9TextileVersion4'
              };

%tex_name = {'G1EpoxyRasterPlate'};

    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
       end
    end
    
    
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
           if exp_tf==1
              rs = exp(100*rs); 
           end  
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0    % no more negative values after the exp()
                rs = rs + abs(texminmin);   % push negative values above zero
           else
                rs = rs - abs(texminmin);   % push lowest value to 0
           end

%            max(rs)
%            min(rs)
           %rs = 0.2 + 0.25.*(1/(texmaxmax-texminmin)).*rs;  % scale input signal to POSFET swing (0.25V), add POSFET offset (0.2V) due to mirror  
           rs = (1/(texmaxmax-texminmin)).*rs;

           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end

    %tsel = 1; % select which texture
    %n = 1; % divide recording into 1second bites; which 1s bite to use
    %L = length(tex_rs{tsel}(n,:));
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),tdelaycnt,mr_cnt,scandir); 
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector       
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    
    
    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling
    

elseif dataset == "LMT_filt"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 108;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate','G1FineAluminumMesh','G1IsolatingFoilMesh','G1MeshFloorCloth', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2','G1ProfiledAluminumPlate','G1RhombAluminumMesh', ...
              'G1RhombAluminumMeshVersion2','G1RoundAluminumMesh','G1RubberMesh','G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta', 'G3AcrylicGlass', 'G3AluminumPlate', ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3EpoxyPlate', 'G3Glass', 'G3StainlessSteel', 'G4Bamboo', ...
              'G4Beech_0', 'G4CherryTree', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4CoreBeech', 'G4Cork', 'G4LaminatedWood', 'G4Larch', ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5RubberPlateVersion2', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2', 'G6IsolatingFoilVersion3', 'G6RedSleaze', ...
              'G6ScouringPadVersion1', 'G6ScouringPadVersion2', 'G6SheepSkin', 'G6Steelwool', ...
              'G6Towel', 'G6VelcroHooks', 'G6VelcroLoops', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7EpdmFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', 'G7MediumFoamVersion2', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8CarbonFoil', 'G8Cardboard', 'G8CardboardVersion2', 'G8GlitterPaperVersion1', ...
              'G8GlitterPaperVersion2', 'G8HighDensityPolyethylen', 'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2', 'G8PlayingCardPaper', 'G8RoughPaper', 'G8Tarp', ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9BlueSatin', 'G9FloorCloth', ...
              'G9FloorClothVersion2', 'G9FloorClothVersion3', 'G9GreenVelvet', 'G9Jeans', ...
              'G9Kashmir', 'G9Leather', 'G9LeatherVersion2', 'G9RedFleeze', ...
              'G9RedVelvet', 'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3', 'G9TextileVersion4'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    
elseif dataset == "LMT_filt_69"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

elseif dataset == "LMT_filt_69_oversampled_10x"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=10;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

elseif dataset == "LMT_filt_69_timediv"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.4;   % prevent MATLAB from crashing + divisible by time_chunks
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/4;  % divide recording to 1 second chunks   
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       trial_new_i=1;
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           %tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(rs);  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [rs,~] = lsim(pos_filt,rs,t);  

           rs_chunk = reshape(rs,length(rs)/4,4);
           rs_chunk = rs_chunk';
           for trial_new=1:4
	   	        tex{k,trial_new_i} = rs_chunk(trial_new,:);
                trial_new_i = trial_new_i + 1;
           end           
       end
    end


       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt*4
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4  
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;   
           %size(rs)
           temp = reshape(rs,[],time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt*4,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

elseif dataset == "LMT_filt_69_oversampled_timediv"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=3.2;   % prevent MATLAB from crashing + divisible by time_chunks
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/4;  % divide recording to 1 second chunks
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       trial_new_i=1;
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           %tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(rs);  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [rs,~] = lsim(pos_filt,rs,t);  

           rs_chunk = reshape(rs,length(rs)/4,4);
           rs_chunk = rs_chunk';
           for trial_new=1:4
	   	        tex{k,trial_new_i} = rs_chunk(trial_new,:);
                trial_new_i = trial_new_i + 1;
           end           
       end
    end


       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt*4
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4  
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;   
           %size(rs)
           temp = reshape(rs,[],time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt*4,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    
elseif dataset == "LMT_filt_69_oversampled_timediv_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 

    
    fs_stim=10e3;
    stim_os=3.2;   % prevent MATLAB from crashing + divisible by time_chunks
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/4;  % divide recording to 1 second chunks
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };

    tex_name = tex_name(K_ind); % get chosen top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       trial_new_i=1;
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           %tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(rs);  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [rs,~] = lsim(pos_filt,rs,t);  

           rs_chunk = reshape(rs,length(rs)/4,4);
           rs_chunk = rs_chunk';
           for trial_new=1:4
	   	        tex{k,trial_new_i} = rs_chunk(trial_new,:);
                trial_new_i = trial_new_i + 1;
           end           
       end
    end


       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt*4
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt*4  
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;   
           %size(rs)
           temp = reshape(rs,[],time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt*4,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

elseif dataset == "LMT_filt_69_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 

    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

elseif dataset == "LMT_filt_69_oversampled_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=2;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    

elseif dataset == "LMT_filt_69_oversampled_10x_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=10;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    
elseif dataset == "LMT_filt_69_oversampled_48x_chunk_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=48;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10; % get only 0.4s out of 4s
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           %size(rs)

           rs = reshape(rs,length(rs)/10,10); % get 1/10th of the sensor signal - first 400ms
           rs = rs';
           tex{k,trial} = rs(1,:);

           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          

           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % factor of 10 to make it right, correct later. time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    


elseif dataset == "LMT_filt_69_oversampled_48x_chunk_topK_gain"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    Av=5;
    fs_stim=10e3;
    stim_os=48;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10; % get only 0.4s out of 4s
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           %size(rs)

           rs = reshape(rs,length(rs)/10,10); % get 1/10th of the sensor signal - first 400ms
           rs = rs';
           tex{k,trial} = rs(1,:);

           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero  
                rs = Av/(texmaxmax-texminmin).*rs;                            
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          

           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % factor of 10 to make it right, correct later. time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling 

elseif dataset == "LMT_69_oversampled_48x_chunk_topK"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=48;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10; % get only 0.4s out of 4s
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = K; %69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           %size(rs)

           rs = reshape(rs,length(rs)/10,10); % get 1/10th of the sensor signal - first 400ms
           rs = rs';
           tex{k,trial} = rs(1,:);

           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           %[tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          

           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % factor of 10 to make it right, correct later. time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling   
    
elseif dataset == "LMT_filt_69_oversampled_48x_chunk"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=48;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10; % get only 0.4s out of 4s
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    %tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           %size(rs)

           rs = reshape(rs,length(rs)/10,10); % get 1/10th of the sensor signal - first 400ms
           rs = rs';
           tex{k,trial} = rs(1,:);

           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          

           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % factor of 10 to make it right, correct later. time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    

elseif dataset == "LMT_69_oversampled_48x_chunk"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=48;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tlength = tlength/10; % get only 0.4s out of 4s
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 69;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1EpoxyRasterPlate', 'G1IsolatingFoilMesh', ...
              'G1PlasticMeshVersion1','G1PlasticMeshVersion2', 'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G1TextileMeshVersion1', 'G2Brick', 'G2CrushedRock', 'G2GraniteTypeVeneziano', ...
              'G2Marble', 'G2RoofTile', 'G2StoneTileVersion1', 'G2StoneTileVersion2', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4Cork', 'G4LaminatedWood',  ...
              'G4ProfiledWoodPlate', 'G4Teak', 'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5FineRubber', 'G5ProfiledRubberPlate', 'G5SolidRubberPlateVersion1', 'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7CoarseFoam', 'G7FineFoamVersion1', 'G7FineFoamVersion2', ...
              'G7FoamFoilVersion1', 'G7FoamPlate', 'G7MediumFoam', ...
              'G7ProfiledFoam', 'G7StyrofoamVersion1', 'G7StyrofoamVersion2', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8PlasticFoilVersion2',  ...
              'G8WallpaperVersion1', 'G8WallpaperVersion2', 'G9FloorCloth', ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', 'G9TextileVersion1', ...
              'G9TextileVersion2', 'G9TextileVersion3'
              };
    
    %tex_name = tex_name(K_ind); % get top K classes

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           rs = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           %size(rs)

           rs = reshape(rs,length(rs)/10,10); % get 1/10th of the sensor signal - first 400ms
           rs = rs';
           tex{k,trial} = rs(1,:);

           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           %[tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          

           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % factor of 10 to make it right, correct later. time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
    

elseif dataset == "LMT_filt_42"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 42;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1IsolatingFoilMesh', ...
              'G1RhombAluminumMesh', ...
              'G1SquaredAluminumMesh', ...
              'G2Brick', ...
              'G2Marble', ...
              'G2StoneTileVersion3', 'G2StoneTypeTerraCotta',  ...
              'G3Brass', 'G3CeramicPlate', 'G3CeramicTile', 'G3Copper', ...
              'G3Glass', 'G3StainlessSteel',  ...
              'G4Beech_0', 'G4CompressedWoodVersion1', 'G4CompressedWoodVersion2', ...
              'G4WoodenPaper', 'G4WoodTypeSilverOak', ....
              'G5SolidRubberPlateVersion2', ...
              'G6Carpet', 'G6CoarseArtificialGrassFibers', 'G6Fibers', 'G6FineArtificialGrassFibers', ...
              'G6IsolatingFoilVersion1', 'G6IsolatingFoilVersion2',  ...
              'G6Steelwool', ...
              'G6Towel', 'G7BumpyFoam', ....
              'G7FineFoamVersion1',  ...
              'G7FoamFoilVersion1', ...
              'G7StyrofoamVersion1', 'G8BubbleFoil', ...
              'G8Cardboard',  ...
              'G8Paper', 'G8PlasticFoilVersion1', ...
              'G8WallpaperVersion2',  ...
              'G9Jeans', ...
              'G9Kashmir', 'G9Leather', ...
              'G9TableClothVersion1', 'G9TableClothVersion2', ...
              'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling  

  elseif dataset == "LMT_filt_11"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=0.5;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 11;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1RhombAluminumMesh', ...
              'G2Brick', ...
              'G2Marble', ...
              'G3Brass', 'G3CeramicPlate', ...
              'G6FineArtificialGrassFibers', ...
              'G7FoamFoilVersion1', ...
              'G8BubbleFoil', ...
              'G8WallpaperVersion2',  ...
              'G9TableClothVersion2', ...
              'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    

   elseif dataset == "LMT_filt_11_oversampled"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=2;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 11;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1RhombAluminumMesh', ...
              'G2Brick', ...
              'G2Marble', ...
              'G3Brass', 'G3CeramicPlate', ...
              'G6FineArtificialGrassFibers', ...
              'G7FoamFoilVersion1', ...
              'G8BubbleFoil', ...
              'G8WallpaperVersion2',  ...
              'G9TableClothVersion2', ...
              'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling    
   elseif dataset == "LMT_filt_11_oversampled_10x"
    %% from LMT dataset - https://dx.doi.org/10.21227/kwsy-x398 
    % used friction data from differential FSR
    %mr_cnt=64;
    
    %fs_stim = 10e3;
    %stim_os = 120;   % to match oversampled Kursun dataset 
    
    fs_stim=10e3;
    stim_os=10;   % prevent MATLAB from crashing
    %DT=1/(fs_stim*stim_os);
    tlength_orig=4.8;   % original time length, based from matrix length NOT FROM PAPER
    scanv_orig = 20; % scanning velocity from paper (mm/s)
    objlength = scanv_orig*tlength_orig;  % length of object scanned
    tlength = objlength/scanv;
    %spatialres = 0.2;   % mm, spatial resolution between sensors
    tdelay = tlength*(spatialres/objlength);  % model space between sensors through time delay
    time_chunks = round(tlength/tdelay);
    tdelaycnt = 2*sqrt(mr_cnt); % no. of delayed samples per texture & trial
    
    disp(['For spatial resolution = ',num2str(spatialres),' mm, scan velocity = ',num2str(scanv),' ,model delay needed is ',num2str(tdelay),' seconds'])
    texcnt = 11;
    %%trialcnt = 10; % 10 recordings per texture
    %trialcnt = 10;
    
    tex_name = {'G1RhombAluminumMesh', ...
              'G2Brick', ...
              'G2Marble', ...
              'G3Brass', 'G3CeramicPlate', ...
              'G6FineArtificialGrassFibers', ...
              'G7FoamFoilVersion1', ...
              'G8BubbleFoil', ...
              'G8WallpaperVersion2',  ...
              'G9TableClothVersion2', ...
              'G9TextileVersion3'
              };

    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    
    for k = 1:numel(tex_name)
       for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr1 = readmatrix([file_loc,'/FSR1/Training/',tex_name{k},'_FrictionFSR1_Movement_train',num2str(trial),'.txt']);
           fsr2 = readmatrix([file_loc,'/FSR2/Training/',tex_name{k},'_FrictionFSR2_Movement_train',num2str(trial),'.txt']);
           tex{k,trial} = fsr1(:,4) - fsr2(:,4);   % differential measurement + ignores weird spaces on dataset file
           
           %% FILTER
           L = length(tex{1,1});  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           [tex{k,trial},~] = lsim(pos_filt,tex{k,trial},t);  
    
           
       end
    end

       
    %% normalize 
%     for k=1:numel(tex_name)
%         for trial=1:trialcnt
%             max(tex{k,trial});
%             texmax(k,trial) = max(tex{k,trial});
%             texmin(k,trial) = min(tex{k,trial});
%         end
%     end
%     texmaxmax = max(max(texmax))
%     texminmin = min(min(texmin))
    
    for k = 1:numel(tex_name)
        for trial=1:trialcnt
           rs = resample(tex{k,trial},stim_os*length(tex{k,trial}),length(tex{k,trial}));   % oversample recording
 
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);                      
        end
    end
    
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))
    
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
           rs = tex{k,trial};
           if exp_tf==0
                rs = rs + abs(texminmin);   % push negative values above zero
                rs = (1/(texmaxmax-texminmin)).*rs;                
           else
                rs = (0.25/(texmaxmax-texminmin)).*rs;      % normalize to expected polymer swing (5N for 50mV/N)              
                rs = 13.4e-6*exp((rs - 0.55)/(1.5*26e-3)); 
           end
           
           tex{k,trial} = rs;
           texmax(k,trial) = max(rs);
           texmin(k,trial) = min(rs);
        end
     end
     
     texmaxmax = max(max(texmax))
     texminmin = min(min(texmin))     
     
     for k = 1:numel(tex_name)
        for trial=1:trialcnt   
     
           rs = tex{k,trial};
            if exp_tf==1   % normalize after exp()
                 rs = (1/(texmaxmax-texminmin)).*rs;
            end     
           
           %tex_rs{k,trial} = rs;          
           temp = reshape(rs,length(rs)/time_chunks,time_chunks); % divide texture sample to groups of tdelay seconds
           temp = temp';
           %tlength_des=4;  % new desired time length
           tlength_des=0.833*tlength;
           for t=1:tdelaycnt     % generate tdelaycnt chunks of texture data (length: tlength_des) offset by tdelay each
                %disp(['--->Creating time chunk # ',num2str(t)]);  
                temp_t = [];
                for i=1:round(tlength_des/tdelay)
                    temp_t = [temp_t,temp(i+t-1,:)];
                end
                tex_rs_t{k,trial,t} = temp_t;   % different receptors get different 1s section (DELAY IN OVERLAP WITH CRF)
           end
       end 
    end
    
%     L = length(tex_rs_t{1}(1,:));  % assume all recording lengths are uniform
%     t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
%     ti = 1:length(t); % time index vector    

    % fill in 64 heminode array, each heminode is offset by tdelay ->
    % simulate up-to-down sweep, assumes pattern repeats on the x-dimension
    %  0  1  2 .... 8
    %  1  2  3 ...  9
    %  2  3  4 ...  10
    %  ...
    %  8  9 10 ...  16
    [tex_rs] = fill_array(tex_rs_t,numel(tex_name),trialcnt,tdelaycnt,mr_cnt,scandir);  
    
    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(scanv_orig/scanv)*(1/(fs_stim*stim_os)); % time vector
    ti = 1:length(t); % time index vector     
    
    %% SPEED UP SIMULATION BY NOT COMPUTING FOR SENSOR DERIVATIVE
%     % differentiate & rectify sensor voltage
%     for k = 1:length(tex_rs)
%         disp(['Reading stimulus # ',num2str(k)]);
%        for h=1:mr_cnt
%             disp(['--->Filling heminode # ',num2str(h),' (diff)']);
%             temp = diff(tex_rs{k}(h,:))*(fs_stim*stim_os);
%             pinds_diff = find(temp > 0);
%             ninds_diff = find(temp < 0);  
% 
%             tex_rsdiffp{k}(h,:) = fill_matrix(temp,pinds_diff,ti);
%             tex_rsdiffn{k}(h,:) = abs(fill_matrix(temp,ninds_diff,ti));   
% 
%             tex_rsdiff{k}(h,:) = tex_rsdiffp{k}(h,:) + tex_rsdiffn{k}(h,:);
%             tex_rsdiff{k}(h,:) = (1/(max(tex_rsdiff{k}(h,:))-min(tex_rsdiff{k}(h,:)))).*tex_rsdiff{k}(h,:);  % scale input signal to fit 1 to 0 swing   
%        end
%     end

    stim.DT = (scanv_orig/scanv)/(fs_stim*stim_os);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength_des;
    stim.fs_stim = fs_stim;
    stim.stim_os = stim_os;    

    %% to be used for adding variability during quantizations
    stim.tdelay = tdelay;
    stim.time_chunks = time_chunks;
    stim.tdelaycnt = tdelaycnt;
    stim.tlength_des = tlength_des;
    stim.tlength_origrs = length(rs); % length after re-sampling     
elseif dataset == "Kylberg_filt_6"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    fs_stim=10e3; % assumption
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = imread([dataset_dir,Files(trial+2).name]);
           rs = rs(1:sqrt(mr_cnt),1:sqrt(mr_cnt));  % get mr_Cnt number of pixels == taxels5
           L = tsamples;  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           t50 = 0.1; % time at 50% rise time
           Tr = 0.1; % risetime
           for h=1:mr_cnt
               disp(['--->Filling receptors # ',num2str(h)]);
               %rs_step(h,:) = [zeros(1,tsamples/4) double(rs(h))*ones(1,tsamples/2) zeros(1,tsamples/4)];  % populate heminodes with step signal with each maximum value equivalent to the pixel value
               rs_step(h,:) =  double(rs(h))./(1+exp(-4*(t-t50)/Tr));               
               %% FILTER

               %[rs_step(h,:),~] = lsim(step_lpf,rs_step(h,:),t); % apply LPF to non-idealize the steps
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_step(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = tsamples; % le

elseif dataset == "Kylberg_filt_6_with_rotations"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    fs_stim=10e3; % assumption
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    
    tex_name = {'canvas1', 'cushion1', ...
              'linseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'-with-rotations/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = imread([dataset_dir,Files(trial+2).name]);
           rs = rs(1:sqrt(mr_cnt),1:sqrt(mr_cnt));  % get mr_Cnt number of pixels == taxels5
           L = tsamples;  % assume all recording lengths are uniform
           t = (0:L-1)*(1/(fs_stim)); % time vector
           t50 = 0.1; % time at 50% rise time
           Tr = 0.1; % risetime
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_step(h,:) = [zeros(1,tsamples/4) double(rs(h))*ones(1,tsamples/2) zeros(1,tsamples/4)];  % populate heminodes with step signal with each maximum value equivalent to the pixel value
               rs_step(h,:) =  double(rs(h))./(1+exp(-4*(t-t50)/Tr));               
               %% FILTER

               %[rs_step(h,:),~] = lsim(step_lpf,rs_step(h,:),t); % apply LPF to non-idealize the steps
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_step(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = tsamples; % le


elseif dataset == "Kylberg_filt_6_scan"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    fs_stim=10e3;
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling
    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_rs(h,:) = resample(rs_scan(h,:),stim_os*length(rs_scan(h,:)),length(rs_scan(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array %tsamples; % le    

elseif dataset == "Kylberg_filt_6_scan_actualdim"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array
elseif dataset == "Kylberg_filt_6_scan"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    fs_stim=10e3;
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling
    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_rs(h,:) = resample(rs_scan(h,:),stim_os*length(rs_scan(h,:)),length(rs_scan(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array %tsamples; % le    

elseif dataset == "Kylberg_filt_6_scan_actualdimscale"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=1;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array

elseif dataset == "Kylberg_filt_6_scan_actualdimscale_chip"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=140;  % gain (small signal due to HPF on step)
    sensor_swingpp = 0.0072;
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               %[x,y] = linto2D(h,mr_cnt);
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               xdim=16; ydim=12;
               [x,y] = linto2D_customxy(h,xdim,ydim);  % reflect actual chip sensors
               %x=x+1; y=y+1;  % x,y starts at 0,0
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan(h,:) = rs(x:x+origx-xdim,y);   % reflect actual chip sensors
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);   

               %tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1               
               tex_rs{k,trial}(h,:) = Av.*(sensor_swingpp/2^8).*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
           end     
       end
    end

    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array

elseif dataset == "Kylberg_filt_6_scan_actualdimscale_chip_lchpf"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 1.956; %0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]);
    %pos_filt = tf([1e3/(2*pi) 0],[1 1007/(2*pi) 7e3/(4*pi^2)]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=140;  % gain (small signal due to HPF on step)
    sensor_swingpp = 0.0072;
    whpf=2*pi*164e-3;
    wlpf=2*pi*43.7*3.5; %2*pi*43.7;  % 3.5 is a tuned number to reflect vout peak2peak
    % closed loop: whpf=10mHz wlpf=1.67kHz
    lc_filt = tf([wlpf 0],[1 whpf+wlpf whpf*wlpf]);
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               %[x,y] = linto2D(h,mr_cnt);
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               xdim=16; ydim=12;
               [x,y] = linto2D_customxy(h,xdim,ydim);  % reflect actual chip sensors
               %x=x+1; y=y+1;  % x,y starts at 0,0
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan(h,:) = rs(x:x+origx-xdim,y);   % reflect actual chip sensors
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t); 
               [rs_step_filt(h,:),~] = lsim(lc_filt,rs_step_filt(h,:),t);

               %tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1               
               tex_rs{k,trial}(h,:) = Av.*(sensor_swingpp/2^8).*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
           end     
       end
    end

    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array

%% replicate VPOS_AFE plot with 7.2mVpp input. Ts=0.1s -> single channel only is activated
elseif dataset == "Kylberg_filt_6_scan_actualdimscale_chip_voutlpf"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.1; %1.956; %0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]); % use this not below since this was used on the code
    %pos_filt = tf([1e3/(2*pi) 0],[1 1007/(2*pi) 7e3/(4*pi^2)]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=140;  % gain (small signal due to HPF on step)
    sensor_swingpp = 0.0072;
    wlpf=2*pi*43.7*3.5;  %84.05Hz (FF); 43.7Hz (TT)
    vout_filt = tf([wlpf],[1 wlpf]);
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               %[x,y] = linto2D(h,mr_cnt);
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               xdim=16; ydim=12;
               [x,y] = linto2D_customxy(h,xdim,ydim);  % reflect actual chip sensors
               %x=x+1; y=y+1;  % x,y starts at 0,0
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan(h,:) = rs(x:x+origx-xdim,y);   % reflect actual chip sensors
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t); 
               [rs_step_filt(h,:),~] = lsim(vout_filt,rs_step_filt(h,:),t);

               %tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1               
               tex_rs{k,trial}(h,:) = Av.*(sensor_swingpp/2^8).*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
           end     
       end
    end

    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array


elseif dataset == "Kylberg_filt_6_scan_actualdimscale180mvpp_chip"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=140;  % gain (small signal due to HPF on step)
    sensor_swingpp = 0.180;
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               %[x,y] = linto2D(h,mr_cnt);
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               xdim=16; ydim=12;
               [x,y] = linto2D_customxy(h,xdim,ydim);  % reflect actual chip sensors
               %x=x+1; y=y+1;  % x,y starts at 0,0
               %rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan(h,:) = rs(x:x+origx-xdim,y);   % reflect actual chip sensors
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);   

               %tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1               
               tex_rs{k,trial}(h,:) = Av.*(sensor_swingpp/2^8).*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
           end     
       end
    end

    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array



elseif dataset == "Kylberg_filt_6_scan_oversampled20x_actualdimscale"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % duration of scanning (assumed)
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=20*round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling

    
    tex_name = {'canvas1', 'cushion1', ...
              'linsseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1e3 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=1;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res             
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/(fs_stim*20)); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array


elseif dataset == "Kylberg_filt_6_with_rotations_scan"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling
    
    tex_name = {'canvas1', 'cushion1', ...
              'linseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'-with-rotations/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_rs(h,:) = resample(rs_scan(h,:),stim_os*length(rs_scan(h,:)),length(rs_scan(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %just make it long enough before chunkification for fill_array 


elseif dataset == "Kylberg_filt_6_with_rotations_scan_actualdim"
    %% from Kylberg dataset - https://kylberg.org/kylberg-texture-dataset-v-1-0/
    % converted PNG image data to a pseudo-texture indendation
    % 576 x 576 pixels -> assume 0.2mm per pixel, total area = 115.2mm
    % (area makes sense when comparing with size of lins seeds)

    texsize=26; %mm, based on the linseed size haha
    fs_stim=10e3; %smaller to reduce memory usage
    texcnt = 6;
    tlength = 0.4; % assumption
    tsamples = tlength*fs_stim;   % number of samples across tlength
    origx=576; origy=576; % no. of pixels on orig texture
    spatialres_orig=texsize/origx;
    stim_os=round(tsamples/origx); % resample sequence of pixel value to have approx, fs_stim sampling
    
    tex_name = {'canvas1', 'cushion1', ...
              'linseeds1','sand1', 'seat2', ...
              'stone1', ...
              };
 
    if randtex==1      
        tex_name = tex_name(randtex_indices); % choose a random subset
    else
        tex_name = tex_name(1:tex_cnt); % choose a subset of the dataset if tex_cnt < 108
    end

    %% SETUP SIGNAL FILTER TO EMULATE POSFET RESPONSE
    %pos_filt = tf([1 0],[1 7]);  % based on Dahiya paper IEEESensors2011
    %step_lpf = tf(3,[1,2,3]); %tf([1],[1 2e3 1e3^2]);  % 2 pole @ 1kHz
    pos_filt = tf([1 0],[1 1007 7e3]);  % based on Dahiya paper IEEESensors2011, add pole at 1kHz to model limited POSFET BW
    Av=100;  % gain (small signal due to HPF on step)
    
    for k = 1:numel(tex_name)
       dataset_dir = [file_loc,'/',tex_name{k},'-with-rotations/']
       Files=dir([dataset_dir]);
       for trial=1:trialcnt
           %Files.name
           %Files(trial).name
           disp(['Reading stimulus # ',num2str(k),': ',tex_name{k},'/',Files(trial+2).name]);
           rs = double(imread([dataset_dir,Files(trial+2).name]));
           for h=1:mr_cnt
               %disp(['--->Filling receptors # ',num2str(h)]);
               %rs_scan(h,:) = rs(h:h+origx-sqrt(mr_cnt));
               [x,y] = linto2D(h,mr_cnt);
               rs_scan(h,:) = rs(x:x+origx-sqrt(mr_cnt),y);
               rs_scan_ds(h,:) = downsample(rs_scan(h,:),round(spatialres/spatialres_orig));  % downsample based on actual spatial res
               rs_scan_rs(h,:) = resample(rs_scan_ds(h,:),stim_os*round(spatialres/spatialres_orig)*length(rs_scan_ds(h,:)),length(rs_scan_ds(h,:)));
               %% FILTER

               L = length(rs_scan_rs(h,:));  % assume all recording lengths are uniform
               t = (0:L-1)*(1/(fs_stim*stim_os)); % time vector
               [rs_step_filt(h,:),~] = lsim(pos_filt,rs_scan_rs(h,:),t);             
               tex_rs{k,trial}(h,:) = Av*rs_step_filt(h,:);  
               %tex_rs{k,trial}(h,:) = rs_step(h,:);  
               tex_rs{k,trial}(h,:) = tex_rs{k,trial}(h,:)/(2^8);  % normalize 8bit pixel value to 1
           end     
       end
    end


    L = length(tex_rs{1,1}(1,:));  % assume all recording lengths are uniform
    t = (0:L-1)*(1/fs_stim); % time vector
    ti = 1:length(t); % time index vector     
    

    stim.DT = 1/(fs_stim);  % needed for tau computation
    stim.tex_rs = tex_rs;
%    stim.tex_rsdiff = tex_rsdiff;
    stim.tex_rsdiff = 0;
    stim.tex_name = tex_name;
    stim.L = L;
    stim.t = t;
    stim.ti = ti;
    stim.tlength = tlength;
    stim.fs_stim = fs_stim;
    stim.stim_os = 1;    

    %% to be used for adding variability during quantizations
    stim.tdelay = 0;  % N/A
    stim.time_chunks = 0; % N/A
    stim.tdelaycnt = 0; % N/A
    stim.tlength_des = tlength;
    stim.tlength_origrs = 2*L; %tsamples - just make it long enough before chunkification for fill_array   

end



end
