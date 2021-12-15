clear; close all;

addpath(genpath('../other_functions'));
addpath('I:\iAtlantic_overflow\VIKING20\D1p2'); % <-  data location

regions = 4; % [1 2 3 4 5 6 7 8 10 11];
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p1_calc_Yn_worker(region)
end
