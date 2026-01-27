% sort spike dataset according to spike times
% to be compatible with Pytorch spike dataloaders

function TD = dataset_sort_timestep(folderPath, num_samples_train, num_samples_test)

if not(isfolder([folderPath,'_sorted']))
    mkdir([folderPath,'_sorted'])
end  

copyfile([folderPath,'/train1K.txt'],[folderPath,'_sorted/','train1K.txt'])
copyfile([folderPath,'/test100.txt'],[folderPath,'_sorted/','test100.txt'])
if isfile([folderPath,'/ave_pop_nSpk.txt'])
    copyfile([folderPath,'/ave_pop_nSpk.txt'],[folderPath,'_sorted/','ave_pop_nSpk.txt'])
end
if isfile([folderPath,'/ave_pop_nSpk_on.txt'])
    copyfile([folderPath,'/ave_pop_nSpk_on.txt'],[folderPath,'_sorted/','ave_pop_nSpk_on.txt'])
end   
if isfile([folderPath,'/ave_pop_nSpk_off.txt'])
    copyfile([folderPath,'/ave_pop_nSpk_off.txt'],[folderPath,'_sorted/','ave_pop_nSpk_off.txt'])
end
if isfile([folderPath,'/ave_SER.txt'])
copyfile([folderPath,'/ave_SER.txt'],[folderPath,'_sorted/','ave_SER.txt'])
end
if isfile([folderPath,'/randtex_indices.txt'])
copyfile([folderPath,'/randtex_indices.txt'],[folderPath,'_sorted/','randtex_indices.txt'])
end

if isfile([folderPath,'/tex_rand-reduced_dset-10-1to10-1to10.txt'])
    copyfile([folderPath,'/tex_rand-reduced_dset-10-1to10-1to10.txt'],[folderPath,'_sorted/','tex_rand-reduced_dset-10-1to10-1to10.txt'])
end
if isfile([folderPath,'/trial_rand-reduced_dset-10-1to10-1to10.txt'])
    copyfile([folderPath,'/trial_rand-reduced_dset-10-1to10-1to10.txt'],[folderPath,'_sorted/','trial_rand-reduced_dset-10-1to10-1to10.txt'])
end
if isfile([folderPath,'/tex_rand-reduced_dset-11-1to10-1to10.txt'])
    copyfile([folderPath,'/tex_rand-reduced_dset-11-1to10-1to10.txt'],[folderPath,'_sorted/','tex_rand-reduced_dset-11-1to10-1to10.txt'])
end
if isfile([folderPath,'/trial_rand-reduced_dset-11-1to10-1to10.txt'])
    copyfile([folderPath,'/trial_rand-reduced_dset-11-1to10-1to10.txt'],[folderPath,'_sorted/','trial_rand-reduced_dset-11-1to10-1to10.txt'])
end
if isfile([folderPath,'/tex_rand-reduced_dset-69-1to10-1to10.txt'])
    copyfile([folderPath,'/tex_rand-reduced_dset-69-1to10-1to10.txt'],[folderPath,'_sorted/','tex_rand-reduced_dset-69-1to10-1to10.txt'])
end
if isfile([folderPath,'/trial_rand-reduced_dset-69-1to10-1to10.txt'])
    copyfile([folderPath,'/trial_rand-reduced_dset-69-1to10-1to10.txt'],[folderPath,'_sorted/','trial_rand-reduced_dset-69-1to10-1to10.txt'])
end

if isfile([folderPath,'/tex_rand-reduced_dset-11-1to40-1to40.txt'])sort
    copyfile([folderPath,'/tex_rand-reduced_dset-11-1to40-1to40.txt'],[folderPath,'_sorted/','tex_rand-reduced_dset-11-1to40-1to40.txt'])
end
if isfile([folderPath,'/trial_rand-reduced_dset-11-1to40-1to40.txt'])
    copyfile([folderPath,'/trial_rand-reduced_dset-11-1to40-1to40.txt'],[folderPath,'_sorted/','trial_rand-reduced_dset-11-1to40-1to40.txt'])
end
if isfile([folderPath,'/tex_rand-reduced_dset-69-1to40-1to40.txt'])
    copyfile([folderPath,'/tex_rand-reduced_dset-69-1to40-1to40.txt'],[folderPath,'_sorted/','tex_rand-reduced_dset-69-1to40-1to40.txt'])
end
if isfile([folderPath,'/trial_rand-reduced_dset-69-1to40-1to40.txt'])
    copyfile([folderPath,'/trial_rand-reduced_dset-69-1to40-1to40.txt'],[folderPath,'_sorted/','trial_rand-reduced_dset-69-1to40-1to40.txt'])
end
if isfile([folderPath,'/th_mismatch.txt'])
    copyfile([folderPath,'/th_mismatch.txt'],[folderPath,'_sorted/','th_mismatch.txt'])
end
if isfile([folderPath,'/th_mismatch_h.txt'])
    copyfile([folderPath,'/th_mismatch_h.txt'],[folderPath,'_sorted/','th_mismatch_h.txt'])
end
if isfile([folderPath,'/th_mismatch_l.txt'])
    copyfile([folderPath,'/th_mismatch_l.txt'],[folderPath,'_sorted/','th_mismatch_l.txt'])
end
if isfile([folderPath,'/LSB.txt'])
    copyfile([folderPath,'/LSB.txt'],[folderPath,'_sorted/','LSB.txt'])
end
if isfile([folderPath,'/tau_refr_mismatch.txt'])
    copyfile([folderPath,'/tau_refr_mismatch.txt'],[folderPath,'_sorted/','tau_refr_mismatch.txt'])
end
if isfile([folderPath,'/trefract.txt'])
    copyfile([folderPath,'/trefract.txt'],[folderPath,'_sorted/','trefract.txt'])
end

for samples = 1:num_samples_train
    filename = [folderPath,'/',num2str(samples),'.bs2']
    TD = Read_Ndataset(filename);

    [TD.ts,I] = sort(TD.ts);
    TD.x = TD.x(I);
    TD.y = TD.y(I);
    TD.p = TD.p(I);   

    Encode_Ndataset([folderPath,'_sorted/',num2str(samples)],TD);

end

for samples = 1:num_samples_test
    filename = [folderPath,'/',num2str(60000+samples),'.bs2']
    TD = Read_Ndataset(filename);

    [TD.ts,I] = sort(TD.ts);
    TD.x = TD.x(I);
    TD.y = TD.y(I);
    TD.p = TD.p(I);    

    Encode_Ndataset([folderPath,'_sorted/',num2str(60000+samples)],TD);

end

end