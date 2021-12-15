clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('I:\iAtlantic_overflow\VIKING20\D1p2'); %('H:/Viking20X'); % <- model data location

regions = [1];
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p2_plot_profiles_worker(region)
end

