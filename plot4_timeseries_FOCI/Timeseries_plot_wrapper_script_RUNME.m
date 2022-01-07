clear; close all;

addpath functions
addpath(genpath('../other_functions'));

in_dir = '/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/INALT10_FOCI/';
expe = {'FOCI1.14-II010','FOCI1.14-II011','FOCI1.14-JH027','FOCI1.19-JH037','FOCI1.19-JH039','FOCI1.14-SW128'};

regions = [1:12]; % 12

for ii= 1 %:length(expe)
    addpath([in_dir expe{ii}])
for aa= 1 %:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p4_plot_timeseries_worker(in_dir,expe{ii},region)
end
end

