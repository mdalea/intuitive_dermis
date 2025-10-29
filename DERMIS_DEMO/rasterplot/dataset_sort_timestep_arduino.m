% sort spike dataset according to spike times
% to be compatible with Pytorch spike dataloaders

%function TD = dataset_sort_timestep_arduino(folderPath, num_samples_train)
function TD = dataset_sort_timestep_arduino(folderPath, samples)

if not(isfolder([folderPath,'_sorted']))
    mkdir([folderPath,'_sorted'])
end  



%for samples = 1:num_samples_train
    filename = [folderPath,'/',num2str(samples),'.bs2']
    %TD = Read_Ndataset(filename);
    TD = Read_Ndataset_longTs(filename);

    [TD.ts,I] = sort(TD.ts);
    TD.x = TD.x(I);
    TD.y = TD.y(I);
    TD.p = TD.p(I);   

    %Encode_Ndataset([folderPath,'_sorted/',num2str(samples)],TD);
    Encode_Ndataset_longTs([folderPath,'_sorted/',num2str(samples)],TD);
%end


end