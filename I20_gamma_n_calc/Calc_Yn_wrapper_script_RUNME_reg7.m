clear; close all;

global_dir = '/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/';
sim_id = 'INALT20_JRA_OMIP/';
in_files = dir([global_dir,sim_id,'monthly/1_INALT20','*reg07_grid_T.nc']);

out_dir = [global_dir,sim_id,'gamma_n/'];
if exist(out_dir,'dir')==0; mkdir(out_dir);end

for idx= 38:length(in_files)
    p1_calc_Yn_worker([global_dir,sim_id,'monthly/'],in_files(idx).name,out_dir);
end
