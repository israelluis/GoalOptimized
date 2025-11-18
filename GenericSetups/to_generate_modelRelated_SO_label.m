%% obtained from running setup quickSO.xml;
folder_path='C:\Users\movea\Dropbox\Collaboration\Guna\GenericSetups';
data_loaded=importdata(fullfile(folder_path,'2_scaled_rajagopal2022_StaticOptimization_force.sto'));
SO_label=data_loaded.colheaders;

save(fullfile(folder_path,'SO_label_rajagopal.mat'),'SO_label');
%% obtained from running setup quickSO.xml;
folder_path   ='C:\Users\movea\Dropbox\Collaboration\Guna\GenericSetups\sample';
folder_outpath='C:\Users\movea\Dropbox\Collaboration\Guna\GenericSetups';
data_loaded=importdata(fullfile(folder_path,'2_scaled_rajagopal2022_StaticOptimization_force.sto'));
SO_label=data_loaded.colheaders;

save(fullfile(folder_outpath,'SO_label_2392.mat'),'SO_label');