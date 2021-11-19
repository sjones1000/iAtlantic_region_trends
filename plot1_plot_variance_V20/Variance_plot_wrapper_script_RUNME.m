clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('H:/Viking20X'); % <- model data location

regions = [2 7 8];
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    disp(['********** Working on Region ' num2str(region) ' ***************']);
    p1_plot_variance_worker(region)
end
