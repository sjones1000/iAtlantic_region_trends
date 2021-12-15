clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/VIKING20X_JRA_OMIP/monthly'); %('H:/Viking20X'); % <- model data location

regions = [1 2 3 4 5 6 7 8 10 11];
for aa= 1%:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p2_plot_TS_worker(region)
end

