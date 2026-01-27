%% generate model spikes for LCS channel
% TD.x - 0 to x (var.)
% TD.y - 0 to y (var.)
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
crf=1;
wta=0;

% -----------------------------
% SWEEPS (edit these arrays)
% -----------------------------
spatialres_array = [0.2];   % example spatial resolutions
mr_cnt_array     = [14*14];   % example MR counts

% >>> NEW SWEEPS <<<
hem_cnt_array    = [4  8 16 64];
global_crf_array = [0 1];
unidist_array    = [0 1];

encode2D=1;
mr_index = [];
hem_coverage = [];
read_mr_index=0;   % 1 if re-use existing heminode branch config, 0 - regenerate random
random_mr=0;
cr=1;   % ASSUMES reusing existing MR indices when testing compression rate. This specifies how much of the mr_index{i} is used per i. Prunes last (1-cr) % of MRs. if cr=1, no compression. TESTED ONLY FOR GLOBAL UNIDIST!!!
%exp_tf=1; % 1 if POSFET converts gate voltage to exp current
tex_cnt_orig=108;
tex_cnt=6;  % if <108, choose only a subset of the dataset for easier classification
trialcnt=10; %40;  % how many trials per class to use
randorder=0;  % 1 if you want order of textures and trials to be random -> make sure to separate and shuffle dataset during pytorch training, 2 if you shuffle list after
randtex=0;  % if 1 - select random tex_cnt from total classes
tex_dataset = 'Kylberg_filt_6_scan_oversampled20x_actualdimscale';
synfilt_en=0;  % 1 - synaptic filter (1st order LPF) on input sensor signal enabled

% make sure base output dirs exist
if ~isfolder([current_loc,'global_outfile'])
    mkdir([current_loc,'global_outfile']);
end
if ~isfolder([current_loc,'global_outfile/spikegen'])
    mkdir([current_loc,'global_outfile/spikegen']);
end
if ~isfolder([current_loc,'global_outfile/spikegen/mr_index'])
    mkdir([current_loc,'global_outfile/spikegen/mr_index']);
end

% -------------------------------------------------------------------------
% IMPORTANT: sweep loops inserted BEFORE mr_cnt/spatialres are used anywhere
% -------------------------------------------------------------------------
for mr_cnt = mr_cnt_array
for spatialres = spatialres_array
for hem_cnt = hem_cnt_array
for global_crf = global_crf_array
for unidist = unidist_array

% reset per sweep (prevents carry-over)
mr_index = [];
hem_coverage = [];

fprintf('SWEEP: mr_cnt=%d, spatialres=%g, hem_cnt=%d, global_crf=%d, unidist=%d\n', ...
    mr_cnt, spatialres, hem_cnt, global_crf, unidist);

%% SKEWED distribution by default
%mr_in_hem used to normalize firing rates
if random_mr==1
    mr_in_hem = ones(1,hem_cnt);
else
    if crf==1
        if global_crf==1 && unidist==1
            mr_in_hem = 8.*ones(1,hem_cnt);   % 8 chosen such that it's sparse enough for global (8/64) and local (8/16) CRFs
            %%mr_in_hem = 16.*ones(1,hem_cnt);   %% COMMENT THIS!!!
        elseif global_crf==1 && unidist == 0
            mr_in_hem = [8 5 3 1];  % factor of 2 to approximate 40MCs/FAI. skewed arbor not really proven yet on MC-FAI. from Lesniak mean: 4-5
            %mr_in_hem = 2.*[11 7 2 0]; % factor of 2 to approximate 40MCs/FAI. skewed arbor not really proven yet on MC-FAI. alternative skewed distribution in Lesniak
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
            hem_coverage(h,:) = [1:mr_cnt];  % choose among all MRs
        end
    else
        i=0;
        for p=[1 5 33 37]  % divide MR array into 4 quadrants
            i=i+1;
            hem_coverage(i,:) = [p:p+3 p+8:p+8+3 p+16:p+16+3 p+24:p+24+3]; % choose among possible MR indices per quadrant
        end
    end

    hem_order = mod(randperm(hem_cnt),4)+1;   % randomize which quadrant get the largest receptor count
    for h=1:hem_cnt
        if global_crf == 1
            idx = randperm(mr_cnt);
        else
            idx = randperm(mr_cnt/4);   % divide into quadrants
        end
        if mod(h,4)==1
            hem_order_2 = mod(randperm(hem_cnt),4)+1;    % randomize which heminode has the receptor count (wraps to 4)
        end
        mr_index{h} = hem_coverage(hem_order(h),idx(1:mr_in_hem(hem_order_2(h))));   % fill in the heminode with the most receptors first, but randomize which quadrant is the heminode
    end
    cnt=hem_cnt;
else
    cnt=mr_cnt;
end

if crf==1 && global_crf==0
   disp(['QI: ']);
   hem_coverage(1,:)
   disp(['QII: ']);
   hem_coverage(2,:)
   disp(['QIII: ']);
   hem_coverage(3,:)
   disp(['QIV: ']);
   hem_coverage(4,:)
end

for h=1:hem_cnt
    mr_index{h}
end

scanv=20; scandir=0;
randtex_indices = randperm(tex_cnt_orig,tex_cnt);   % choose a random subset

[stim] = read_textures('../../tactile_dataset/Kylberg-6',tex_dataset,tex_cnt,trialcnt,randtex,randtex_indices,fs_stim,stim_os,mr_cnt,scanv,scandir,0,spatialres);

% ensure base output dir exists (again, in case path changes)
if ~isfolder([current_loc,'global_outfile/spikegen'])
    mkdir([current_loc,'global_outfile/spikegen']);
end

tex_name = stim.tex_name;
writecell(tex_name,[current_loc,'global_outfile/spikegen/tex_name.txt']);
[maxSR,sig_ave,sig_pp] = print_textures(current_loc,stim,'SA',plot_trials,mr_cnt);

mean(sig_pp) %mean peak to peak per class
mean(mean(sig_pp)) %mean peak to peak of signal (from 0.1s to 2s)

noiseVar = [0 0 0 0 0];

noiseVar_dim = size(noiseVar(:,:));

trials_train=round(tex_cnt*trialcnt*0.7);  trials_test=round(tex_cnt*trialcnt*0.3);
trials_for_train=[1:trialcnt];  trials_for_test=[1:trialcnt];
trials = trials_train + trials_test;

% Generate order of texture presentation
if randorder == 0 || randorder==2
    tex_rand = [];
    for texture=1:tex_cnt
        tex_rand = [tex_rand texture*ones(1,trialcnt)];
    end
else
    tex_rand = randi(numel(tex_name),1,trials);  % generate random order of texture presentations as training data
end

% size changed from cnt to mr_cnt to preserve orig size; sums happen on
% LCSDQ later after adding noise
for i=1:noiseVar_dim(1)
    th_mismatch_h = sqrt(noiseVar(i,4))*randn(1,mr_cnt); % generate mismatch spread on all taxel thresholds
    th_mismatch_l = sqrt(noiseVar(i,4))*randn(1,mr_cnt); % generate mismatch spread on all taxel thresholds

    for N_th=N_th_array
        if crf == 1
            LSB_th = max(mr_in_hem)*(1/(2^N_th));  % higher threshold to normalize firing rate
            th_mismatch = LSB_th*noiseVar(i,4)*randn(1,mr_cnt);
        else
            LSB_th = 1/(2^N_th);
            th_mismatch = LSB_th*noiseVar(i,4)*randn(1,mr_cnt);
        end

        for tt=1:size(trials_for_train,1)
            if randorder == 0 || randorder == 2
                trial_rand = [];
                for texture=1:tex_cnt
                    trial_rand = [trial_rand 1:1:trialcnt];
                end
                if randorder == 2
                   iidx = randperm(length(tex_rand));
                   % shuffle texture and trial order together
                   tex_rand = tex_rand(iidx);
                   trial_rand = trial_rand(iidx);
                end
            else
                trial_rand_train = randi([min(trials_for_train(tt,:)) max(trials_for_train(tt,:))],1,trials_train);
                trial_rand_test = randi([min(trials_for_test(tt,:)) max(trials_for_test(tt,:))],1,trials_test);
                trial_rand = [trial_rand_train trial_rand_test];
            end

            for exp_tf = 0

                folderPath = [current_loc,'global_outfile/spikegen/ALL_N_lcadc-',num2str(1e6*noiseVar(i,1)),'-uV2-',num2str(1e6*noiseVar(i,2)),'-uV2-',num2str(1e6*noiseVar(i,3)),'-uV2-',num2str(noiseVar(i,5)),'-mult-',num2str(1e6*noiseVar(i,4)),'-uV2-',num2str(N_th)];
                folderPath_mrindex = [current_loc,'global_outfile/spikegen/mr_index'];

                folderPath_ext = [];

                if random_mr~=1
                    if crf==1
                        folderPath_ext = [folderPath_ext,'-crf-',num2str(hem_cnt)];
                    end

                    % make global/local explicit, but preserve old naming for global=1
                    if global_crf==1
                        folderPath_ext = [folderPath_ext,'-global_crf'];
                    else
                        folderPath_ext = [folderPath_ext,'-local_crf'];
                    end

                    % make unidist/skew explicit, but preserve old naming for unidist=1
                    if unidist==1
                        folderPath_ext = [folderPath_ext,'-unidist-max',num2str(max(mr_in_hem))];
                    else
                        folderPath_ext = [folderPath_ext,'-skewdist'];
                    end

                    % NOTE: mr_cnt tag removed from folderPath_ext so it does NOT affect folderPath_mrindex
                    % (mr_cnt will be appended ONLY to folderPath later)
                else
                    folderPath_ext = [folderPath_ext,'-random_mr-',num2str(hem_cnt)];  % pretend hem_cnt is MR_cnt
                end

                folderPath = [folderPath,folderPath_ext];
                folderPath_mrindex = [folderPath_mrindex,folderPath_ext];

                if wta==1
                    folderPath = [folderPath,'-wta'];
                end

                if cr~=1
                    folderPath = [folderPath,'-cr-',num2str(100*cr),'pct'];
                    folderPath_mrindex = [folderPath_mrindex,'-cr-',num2str(100*cr),'pct'];
                end

                if exp_tf==1
                    folderPath = [folderPath,'-exp_tf'];
                    folderPath_mrindex = [folderPath_mrindex,'-exp_tf'];
                end

                if size(stim.tex_rs,2) > 1
                    folderPath = [folderPath,'-multitrial'];
                    folderPath_mrindex = [folderPath_mrindex,'-multitrial'];
                end

                if randtex == 1
                    folderPath = [folderPath,'-randomtexture6'];
                    folderPath_mrindex = [folderPath_mrindex,'-randomtexture6'];
                end

                if (tex_dataset == "LMT_filt" || tex_dataset == "LMT") && tex_cnt < 108
                    folderPath = [folderPath,'-reduced_dset-',num2str(tex_cnt)];
                    folderPath_mrindex = [folderPath_mrindex,'-reduced_dset-',num2str(tex_cnt)];
                end

                if (tex_dataset == "LMT_filt_69" || tex_dataset == "LMT_69" || tex_dataset == "LMT_filt_69_oversampled" || tex_dataset == "LMT_filt_69_oversampled_10x" || tex_dataset == "LMT_filt_69_timediv" || tex_dataset == "LMT_filt_69_oversampled_48x_chunk") && tex_cnt < 69
                    folderPath = [folderPath,'-reduced_dset-',num2str(tex_cnt)];
                    folderPath_mrindex = [folderPath_mrindex,'-reduced_dset-',num2str(tex_cnt)];
                end

                if (tex_dataset == "LMT_filt_11" || tex_dataset == "LMT_11" || tex_dataset == "LMT_filt_11_oversampled" || tex_dataset == "LMT_filt_11_oversampled_10x") && tex_cnt < 11
                    folderPath = [folderPath,'-reduced_dset-',num2str(tex_cnt)];
                    folderPath_mrindex = [folderPath_mrindex,'-reduced_dset-',num2str(tex_cnt)];
                end

                if (tex_dataset == "LMT_filt_69_oversampled_timediv_topK" || tex_dataset == "LMT_filt_69_topK" || tex_dataset == "LMT_filt_69_oversampled_topK" || tex_dataset == "LMT_filt_69_oversampled_10x_topK" || tex_dataset == "LMT_filt_69_oversampled_48x_chunk_topK") && tex_cnt < 10
                    folderPath = [folderPath,'-reduced_dset-',num2str(tex_cnt)];
                    folderPath_mrindex = [folderPath_mrindex,'-reduced_dset-',num2str(tex_cnt)];
                end

                if (tex_dataset == "Kursun_filt" || tex_dataset == "Kursun") && tex_cnt < 12
                    folderPath = [folderPath,'-reduced_dset-',num2str(tex_cnt)];
                    folderPath_mrindex = [folderPath_mrindex,'-reduced_dset-',num2str(tex_cnt)];
                end

                folderPath = [folderPath,'-',tex_dataset,'-'];

                if size(stim.tex_rs,2) > 1
                    folderPath = [folderPath,'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2)))];
                    folderPath_mrindex = [folderPath_mrindex,'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2)))];
                end

                if randorder == 1
                    folderPath = [folderPath,'-randomorder'];
                    folderPath_mrindex = [folderPath_mrindex,'-randomorder'];
                end

                if randorder == 2
                    folderPath = [folderPath,'-randomorder2'];
                    folderPath_mrindex = [folderPath_mrindex,'-randomorder2'];
                end

                % ---------------------------------------------------------
                % Append sweep params ONLY to folderPath (NOT folderPath_mrindex)
                % ---------------------------------------------------------
                folderPath = [folderPath,'-mr_cnt-',num2str(mr_cnt),'-spatialres-',num2str(spatialres)];

                if not(isfolder(folderPath))
                    mkdir(folderPath);
                end
                if not(isfolder(folderPath_mrindex))
                    mkdir(folderPath_mrindex);
                end

                % moved here (after folder creation)
                if crf==1
                    if read_mr_index == 1 || cr~=1   % always assume reusing MR indices when testing spatial compression rate
                        orig_mr_index_numel = numel(mr_index);
                        for j=1:orig_mr_index_numel
                            mr_index{j} = readmatrix([folderPath_mrindex,'/mr_index',num2str(j),'.txt']);
                            mr_index{j} = mr_index{j}(1:round(cr*numel(mr_index{j})));
                        end
                    else
                        for j=1:numel(mr_index)
                            mr_index{j} = mr_index{j}(1:round(cr*numel(mr_index{j})));

                            % write both (central cache + per-run folder), so reuse works and local debugging is easy
                            writematrix(mr_index{j},[folderPath_mrindex,'/mr_index',num2str(j),'.txt']);
                            writematrix(mr_index{j},[folderPath,'/mr_index',num2str(j),'.txt']);
                        end
                    end
                end

                fid_tr = fopen([folderPath,'/train1K.txt'],'wt');
                fid_te = fopen([folderPath,'/test100.txt'],'wt');

                fprintf(fid_tr,'#sample #class\n');
                fprintf(fid_te,'#sample #class\n');

                tot_nSpk_on=0;   tot_nSpk_off=0;
                tot_SER=0;

                mr_in_hem = round(cr.*mr_in_hem);  % actual number of MR in hem truncated. For the quantization code to see
                max_maxSR = max(max(max(maxSR)));

                if retrieve==1
                    th_mismatch = readmatrix([folderPath,'/th_mismatch.txt']);
                    randtex_indices = readmatrix([folderPath,'/randtex_indices.txt']);
                else
                    writematrix(th_mismatch,[folderPath,'/th_mismatch.txt']);
                    writematrix(randtex_indices,[folderPath,'/randtex_indices.txt']);
                end

                % save texture & trial order just in case Matlab crashes so we can restart
                if retrieve==1
                    if numel(tex_name) < 108
                        trial_rand = readmatrix([folderPath,'/trial_rand-reduced_dset-',num2str(numel(tex_name)),'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                        tex_rand = readmatrix([folderPath,'/tex_rand-reduced_dset-',num2str(numel(tex_name)),'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                    else
                        trial_rand = readmatrix([folderPath,'/trial_rand-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                        tex_rand = readmatrix([folderPath,'/tex_rand-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                    end
                else
                    if numel(tex_name) < 108
                        writematrix(trial_rand,[folderPath,'/trial_rand-reduced_dset-',num2str(numel(tex_name)),'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                        writematrix(tex_rand,[folderPath,'/tex_rand-reduced_dset-',num2str(numel(tex_name)),'-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                    else
                        writematrix(trial_rand,[folderPath,'/trial_rand-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                        writematrix(tex_rand,[folderPath,'/tex_rand-',num2str(trials_for_train(tt,1)),'to',num2str(trials_for_train(tt,size(trials_for_train,2))),'-',num2str(trials_for_test(tt,1)),'to',num2str(trials_for_test(tt,size(trials_for_test,2))),'.txt']);
                    end
                end

                if retrieve==1
                    starti = 161; %1932+618;
                else
                    starti = 1;
                end

                %% if file writing on train.txt & test.txt FAILS, no file is written since it is not closed properly
                %% DO THIS EVERYTIME INSTEAD
                for textures=1:numel(tex_rand)
                     if textures <= trials_train
                        if encode2D == 0
                            disp(['Writing to ',folderPath,'/',num2str(textures),'.bs1 : stimuli = ',tex_name{tex_rand(textures)}]);
                        elseif encode2D == 1
                            disp(['Writing to ',folderPath,'/',num2str(textures),'.bs2 : stimuli = ',tex_name{tex_rand(textures)}]);
                        end
                        fprintf(fid_tr,'%d\t%d\n',textures,tex_rand(textures)-1);    % append
                    else
                         if encode2D == 0
                            disp(['Writing to ',folderPath,'/',num2str(60000+textures-trials_train),'.bs1 : stimuli = ',tex_name{tex_rand(textures)}]);
                        elseif encode2D == 1
                            disp(['Writing to ',folderPath,'/',num2str(60000+textures-trials_train),'.bs2 : stimuli = ',tex_name{tex_rand(textures)}]);
                        end
                        fprintf(fid_te,'%d\t%d\n',60000+textures-trials_train,tex_rand(textures)-1);       % append
                    end
                end

                for textures = starti:numel(tex_rand)
                    tloop_lc = (1/(2^N_th)) / max_maxSR;
                    hem = [];

                    hem = LC_Quantization(current_loc,tex_rand(textures),trial_rand(textures),N_th,stim,max_maxSR,'SA',mr_index,mr_in_hem,crf,global_crf,wta,mr_cnt,hem_cnt,noiseVar(i,:),th_mismatch_h,th_mismatch_l,scandir,fileprint);  % contains SER,sampcnt,maxSR,tlooprequired,and spktimes

                    % clear
                    eventStamp = [];
                    TD.x = []; TD.y = []; TD.p = []; TD.ts = [];

                    if wta==1    % new cnt for printing
                        cnt = 1;   % 1 node for now
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

                        if nSpk_on == 0   % handle case for empty spike train since SNN on Python does not like empty tensors
                            if encode2D == 0
                                newrow = [h   0    1];
                                eventStamp = [eventStamp;newrow];
                            elseif encode2D == 1
                                [x,y] = linto2D_customxy(h,16,12);
                                x= x-1; y=y-1;
                                TD.x = [TD.x;x];  TD.y = [TD.y;y];
                                TD.p = [TD.p;1];
                                TD.ts = [TD.ts;0];
                            end
                        else
                            for s=1:nSpk_on
                                if encode2D == 0
                                    newrow = [h   tempspike_on(s)*1e6    1];
                                    eventStamp = [eventStamp;newrow];
                                elseif encode2D == 1
                                    [x,y] = linto2D_customxy(h,16,12);
                                    x= x-1; y=y-1;
                                    TD.x = [TD.x;x];  TD.y = [TD.y;y];
                                    TD.p = [TD.p;1];
                                    TD.ts = [TD.ts  ;  cast( tempspike_on(s)*1e6,  'uint32') ];
                                end
                            end
                        end

                        % OFF spikes
                        nSpk_off = numel(hem.spktimes_off{h});
                        tot_nSpk_off = tot_nSpk_off + nSpk_off;
                        tempspike_off = hem.spktimes_off{h};
                        tempspike_off(tempspike_off < 0) = 0;

                        if nSpk_off == 0
                            if encode2D == 0
                                newrow = [h   0    1];
                                eventStamp = [eventStamp;newrow];
                            elseif encode2D == 1
                                [x,y] = linto2D_customxy(h,16,12);
                                x= x-1; y=y-1;
                                TD.x = [TD.x;x];  TD.y = [TD.y;y];
                                TD.p = [TD.p;2];
                                TD.ts = [TD.ts;0];
                            end
                        else
                            for s=1:nSpk_off
                                if encode2D == 0
                                    newrow = [h   tempspike_off(s)*1e6    1];
                                    eventStamp = [eventStamp;newrow];
                                elseif encode2D == 1
                                    [x,y] = linto2D_customxy(h,16,12);
                                    x= x-1; y=y-1;
                                    TD.x = [TD.x;x];  TD.y = [TD.y;y];
                                    TD.p = [TD.p;2];
                                    TD.ts = [TD.ts  ;  cast( tempspike_off(s)*1e6,  'uint32') ];
                                end
                            end
                        end
                    end

                    if textures <= trials_train
                        if encode2D == 0
                            disp(['Writing to ',folderPath,'/',num2str(textures),'.bs1 : stimuli = ',tex_name{tex_rand(textures)}]);
                            encode1DBinSpikes([folderPath,'/',num2str(textures),'.bs1'],eventStamp);
                        elseif encode2D == 1
                            disp(['Writing to ',folderPath,'/',num2str(textures),'.bs2 : stimuli = ',tex_name{tex_rand(textures)}]);
                            Encode_Ndataset([folderPath,'/',num2str(textures)],TD);
                        end
                    else
                         if encode2D == 0
                            disp(['Writing to ',folderPath,'/',num2str(60000+textures-trials_train),'.bs1 : stimuli = ',tex_name{tex_rand(textures)}]);
                            encode1DBinSpikes([folderPath,'/',num2str(60000+textures-trials_train),'.bs1'],eventStamp);
                        elseif encode2D == 1
                            disp(['Writing to ',folderPath,'/',num2str(60000+textures-trials_train),'.bs2 : stimuli = ',tex_name{tex_rand(textures)}]);
                            Encode_Ndataset([folderPath,'/',num2str(60000+textures-trials_train)],TD);
                        end
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

                dataset_sort_timestep(folderPath,trials_train,trials_test);  % sort spikes for binning in pytorch

            end % exp_tf
        end % tt
    end % N_th
end % noiseVar

end % unidist sweep
end % global_crf sweep
end % hem_cnt sweep
end % spatialres sweep
end % mr_cnt sweep

toc
toc
