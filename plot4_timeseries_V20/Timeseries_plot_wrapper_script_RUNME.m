clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('I:\iAtlantic_overflow\VIKING20\D1p2'); %  addpath('H:/Viking20X'); % <-  data location

regions = 2; % [1:12]; % 12
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(region)
end

