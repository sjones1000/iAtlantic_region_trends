clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('H:/EN4/iAtlantic'); % <-  data location

regions = [1 12]; %[1:12]; % 12
for aa= 1:length(regions)
    region = regions(aa); 
    %p2_plot_TS_worker(region)
    p2_plot_TS_worker_colour_by_pres(region)
end

