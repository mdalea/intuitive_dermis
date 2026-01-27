%% Load texture dataset and print

tic

clear all
clearvars

%current_loc = '/users/micas/malea/Documents/all_scripts-slayer/';
%cd '/users/micas/malea/Documents/all_scripts-slayer/'

fs_stim = 200; % original signal frequency -> will be overwritten in read_textures
tdelay = 100e-9; % intended loop delay (how often crossing can be detected)
stim_os = cast(1/(tdelay*fs_stim),'double');
stim_os = 6000;

mr_cnt=12*16;
%tex_dataset = 'Kylberg_filt_6_scan_actualdimscale_chip_lchpf'
tex_dataset = 'Kylberg_filt_6_scan_actualdimscale_chip_voutlpf'
tex_cnt=6;  % classes
trialcnt=1 %40;  % how many trials per class to use
randorder=0;  % 1 if you want order of textures and trials to be random -> make sure to separate and shuffle dataset during pytorch training
randtex=0;  % if 1 - select random tex_cnt from total classes
randtex_indices = randperm(tex_cnt,tex_cnt)   % choose a random subset  (NOT USED HERE)
synfilt_en=0;  % 1 - synaptic filter (1st order LPF) on input sensor signal enabled
spatialres=0.2; % mm, spatial resolution of the sensors (needed by read_texture)

scanv=20; scandir=0;
[stim_orig] = read_textures_arduino([pwd,'/tactile_dataset/Kylberg-6'],tex_dataset,tex_cnt,trialcnt,randtex,randtex_indices,fs_stim,stim_os,mr_cnt,scanv,scandir,0,spatialres);     
[maxSR,sig_ave] = print_textures(pwd,stim_orig,'SA',0,mr_cnt);