%% generate model spikes for LCS channel
% TD.x - 0 to 15
% TD.y - 0 to 11
% TD.p - 1 (ON) & 2 (OFF)
% TD.ts - microseconds

tic

clear all
clearvars

current_loc = [pwd,'/'];

fileprint=0;
plot_trials=1; % 1 - print a subset of the trials per class; 0 - print all
retrieve=0;

fs_stim = 200; % original signal frequency
tdelay = 100e-9; % intended loop delay (how often crossing can be detected)
stim_os = cast(1/(tdelay*fs_stim),'double');
stim_os = 6000;

% use as default
N_th_array = [3]; %[2 3 4 5 6 7];
crf=0;
unidist=0;
global_crf=0;
wta=0;
hem_cnt=0;  % ignore this if crf==0
encode2D=1;
mr_index = [];
hem_coverage = [];
read_mr_index=0;   % 1 if re-use existing heminode branch config, 0 - regenerate random
random_mr=0;   
cr=1;   % compression rate
tex_cnt_orig=108;
tex_cnt=6;  % subset of dataset
trialcnt=40; %40;  % how many trials per class
randorder=0;  
randtex=0;  
%tex_dataset = 'Kylberg_filt_6_scan_actualdimscale_chip_lchpf';
tex_dataset = 'Kylberg_filt_6_scan_oversampled20x_actualdimscale'
synfilt_en=0;  % 1 - synaptic filter enabled

%% Sweep parameters
spatialres_array = [0.025 0.05 0.1 0.2 0.4 1 2 5];   % example spatial resolutions
%spatialres_array = [2 5];   % example spatial resolutions
mr_cnt_array = [3*3 8*8 14*14 16*16 18*18];      % example MR counts

for spatialres = spatialres_array
    for mr_cnt = mr_cnt_array

        %% SKEWED distribution by default
        if random_mr==1
            mr_in_hem = ones(1,hem_cnt);
        else
            if crf==1
                if global_crf==1 && unidist==1 
                    mr_in_hem = 8.*ones(1,hem_cnt);  
                elseif global_crf==1 && unidist == 0    
                    mr_in_hem = [8 5 3 1];  
                elseif global_crf==0 && unidist == 1  
                    mr_in_hem = 8.*ones(1,hem_cnt);  
                elseif global_crf==0 && unidist == 0  
                    mr_in_hem = [8 5 3 1];     
                end
            else
               mr_in_hem = 1; 
            end
        end

        %% GENERATE RECEPTOR INDICES GROUPED BY HEMINODE
        if crf==1
            if global_crf==1
                for h=1:hem_cnt
                    hem_coverage(h,:) = [1:mr_cnt];  
                end
            else
                i=0;
                for p=[1 5 33 37]  % divide MR array into 4 quadrants
                    i=i+1;        
                    hem_coverage(i,:) = [p:p+3 p+8:p+8+3 p+16:p+16+3 p+24:p+24+3];
                end
            end

            hem_order = mod(randperm(hem_cnt),4)+1;  
            for h=1:hem_cnt
                if global_crf == 1
                    idx = randperm(mr_cnt);
                else
                    idx = randperm(mr_cnt/4);  
                end
                if mod(h,4)==1
                    hem_order_2 = mod(randperm(hem_cnt),4)+1;   
                end    
                mr_index{h} = hem_coverage(hem_order(h),idx(1:mr_in_hem(hem_order_2(h))));  
            end
            cnt=hem_cnt;
        else
            cnt=mr_cnt;
        end

        %% Read textures
        scanv=20; scandir=0;
        randtex_indices = randperm(tex_cnt_orig,tex_cnt);  

        [stim] = read_textures('../../tactile_dataset/Kylberg-6', ...
            tex_dataset, tex_cnt, trialcnt, randtex, randtex_indices, ...
            fs_stim, stim_os, mr_cnt, scanv, scandir, 0, spatialres);    

        mkdir 'global_outfile/spikegen'
        tex_name = stim.tex_name;  
        writecell(tex_name,[current_loc,'global_outfile/spikegen/tex_name.txt']);    

        [maxSR,sig_ave,sig_pp] = print_textures(current_loc,stim,'SA',plot_trials,mr_cnt);

        mean(sig_pp) 
        mean(mean(sig_pp)) 

        noiseVar = [0 0 0 0 0];
        noiseVar_dim = size(noiseVar(:,:));

        trials_train=round(tex_cnt*trialcnt*0.7);  
        trials_test=round(tex_cnt*trialcnt*0.3);
        trials_for_train=[1:trialcnt];  
        trials_for_test=[1:trialcnt];
        trials = trials_train + trials_test;

        if randorder == 0 || randorder==2
            tex_rand = [];
            for texture=1:tex_cnt
                tex_rand = [tex_rand texture*ones(1,trialcnt)];
            end   
        else
            tex_rand = randi(numel(tex_name),1,trials);
        end

        %% Generate spikes for all noise and thresholds
        for i=1:noiseVar_dim(1)
            th_mismatch_h = sqrt(noiseVar(i,4))*randn(1,cnt); 
            th_mismatch_l = sqrt(noiseVar(i,4))*randn(1,cnt); 

            for N_th=N_th_array
                if crf == 1
                    LSB_th = max(mr_in_hem)*(1/(2^N_th));  
                    th_mismatch = LSB_th*noiseVar(i,4)*randn(1,cnt);
                else
                    LSB_th = 1/(2^N_th);
                    th_mismatch = LSB_th*noiseVar(i,4)*randn(1,cnt);
                end  

                for tt=1:size(trials_for_train,1)
                    trial_rand = 1:trialcnt;
                    exp_tf = 0;

                    folderPath = [current_loc,'global_outfile/spikegen/ALL_N_lcadc-',...
                        num2str(1e6*noiseVar(i,1)),'-uV2-',...
                        num2str(1e6*noiseVar(i,2)),'-uV2-',...
                        num2str(1e6*noiseVar(i,3)),'-uV2-',...
                        num2str(noiseVar(i,5)),'-mult-',...
                        num2str(1e6*noiseVar(i,4)),'-uV2-',...
                        num2str(N_th),...
                        '-mr_cnt-', num2str(mr_cnt),...
                        '-spatialres-', strrep(num2str(spatialres),'.','p')];  % replace '.' with 'p' to make valid folder name

                    folderPath_mrindex = [current_loc,'global_outfile/spikegen/mr_index'];

                    if crf==1 && cr~=1
                        for j=1:numel(mr_index)
                            mr_index{j} = mr_index{j}(1:round(cr*numel(mr_index{j})));
                            writematrix(mr_index{j},[folderPath,'/mr_index',num2str(j),'.txt']);
                        end                
                    end

                    if not(isfolder(folderPath))
                        mkdir(folderPath)
                    end  

                    fid_tr = fopen([folderPath,'/train1K.txt'],'wt');
                    fid_te = fopen([folderPath,'/test100.txt'],'wt');  
                    fprintf(fid_tr,'#sample #class\n');
                    fprintf(fid_te,'#sample #class\n');         

                    tot_nSpk_on=0;   tot_nSpk_off=0;
                    tot_SER=0;

                    max_maxSR = max(max(max(maxSR)));
                    starti = 1; 

                    %% Main spike generation loop
                    for textures = starti:numel(tex_rand)    
                        hem = LC_Quantization(current_loc,tex_rand(textures),trial_rand(textures),...
                            N_th,stim,max_maxSR,'SA',mr_index,mr_in_hem,crf,global_crf,wta,...
                            mr_cnt,hem_cnt,noiseVar(i,:),th_mismatch_h,th_mismatch_l,scandir,fileprint);

                        eventStamp = [];             
                        TD.x = []; TD.y = []; TD.p = []; TD.ts = [];

                        if wta==1
                            cnt = 1;  
                        else
                            if crf==1
                                cnt = hem_cnt;
                            else
                                cnt = mr_cnt;
                            end
                        end

                        for h=1:cnt
                            % ON spikes
                            nSpk_on = numel(hem.spktimes_on{h});
                            tot_nSpk_on = tot_nSpk_on + nSpk_on;
                            tot_SER = tot_SER + hem.SER(h);
                            tempspike_on = hem.spktimes_on{h}; 
                            tempspike_on(tempspike_on < 0) = 0;                      

                            % OFF spikes
                            nSpk_off = numel(hem.spktimes_off{h});
                            tot_nSpk_off = tot_nSpk_off + nSpk_off;
                            tempspike_off = hem.spktimes_off{h}; 
                            tempspike_off(tempspike_off < 0) = 0;                      

                            % Encode spikes into TD or eventStamp
                            for s=1:nSpk_on
                                [x,y] = linto2D_customxy(h,16,12);
                                x= x-1; y=y-1;
                                TD.x = [TD.x;x]; TD.y = [TD.y;y];
                                TD.p = [TD.p;1];
                                TD.ts = [TD.ts; cast(tempspike_on(s)*1e6,'uint32')];
                            end

                            for s=1:nSpk_off
                                [x,y] = linto2D_customxy(h,16,12);
                                x= x-1; y=y-1;
                                TD.x = [TD.x;x]; TD.y = [TD.y;y];
                                TD.p = [TD.p;2];
                                TD.ts = [TD.ts; cast(tempspike_off(s)*1e6,'uint32')];
                            end
                        end

                        if textures <= trials_train        
                            Encode_Ndataset([folderPath,'/',num2str(textures)],TD);        
                        else
                            Encode_Ndataset([folderPath,'/',num2str(60000+textures-trials_train)],TD);        
                        end   
                    end

                    fclose(fid_tr);
                    fclose(fid_te);

                    ave_pop_nSpk_on = tot_nSpk_on/(trials*stim.tlength);
                    ave_pop_nSpk_off = tot_nSpk_off/(trials*stim.tlength);
                    ave_pop_nSpk = ave_pop_nSpk_on + ave_pop_nSpk_off;
                    ave_SER = tot_SER/(cnt*trials);
                    writematrix(ave_pop_nSpk_on,[folderPath,'/ave_pop_nSpk_on.txt']);
                    writematrix(ave_pop_nSpk_off,[folderPath,'/ave_pop_nSpk_off.txt']);
                    writematrix(ave_pop_nSpk,[folderPath,'/ave_pop_nSpk.txt']);
                    writematrix(ave_SER,[folderPath,'/ave_SER.txt']);

                    dataset_sort_timestep(folderPath,trials_train,trials_test);
                end
            end
        end
    end
end

toc

