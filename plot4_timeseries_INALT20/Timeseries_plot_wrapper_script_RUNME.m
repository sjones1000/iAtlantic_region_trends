clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('I:\iAtlantic_overflow\INALT20\INALT20_JRA_OMIP\INALT20_JRA_OMIP'); % <-  data location

regions = [7 8]; % [7 8 9 10 11 12]; 
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(region)
end

