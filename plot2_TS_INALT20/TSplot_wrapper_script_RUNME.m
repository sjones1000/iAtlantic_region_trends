clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('I:\iAtlantic_overflow\INALT20\INALT20_JRA_OMIP\INALT20_JRA_OMIP');  % <- model data location

regions = [7 8 9 10 11 12];
for aa= 1:length(regions)
    region = regions(aa);
    p2_plot_TS_worker(region)
end
