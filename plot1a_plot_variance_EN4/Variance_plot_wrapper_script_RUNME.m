clear; close all;

addpath functions
addpath(genpath('../other_functions'));
addpath('H:/EN4'); % <-  data location

regions =  [1:12];
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    disp(['********** Working on Region ' num2str(region) ' ***************']);
    p1_plot_variance_worker(region)
end
