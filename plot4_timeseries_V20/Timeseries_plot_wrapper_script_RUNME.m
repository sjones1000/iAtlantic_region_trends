clear; close all;

addpath functions
addpath(genpath('../other_functions'));

in_dir = '/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/VIKING20X_JRA_OMIP/';

regions = [1,2,4,5] %[3,6,7,8,10,11]; %4; 
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(in_dir,region)
end

