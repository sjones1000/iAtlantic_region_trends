clear; close all;

addpath functions
addpath(genpath('../other_functions'));

in_dir = '/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/INALT10_FOCI/';

regions = [1:12]; % 12
    
for aa= 4 %2:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(in_dir,region)
end


