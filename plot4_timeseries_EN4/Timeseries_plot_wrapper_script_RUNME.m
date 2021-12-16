clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('H:/EN4/iAtlantic'); % <-  data location

regions = [1:12]; % 12
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(region)
end

