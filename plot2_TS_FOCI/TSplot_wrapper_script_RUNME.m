clear; close all;

in_dir = '/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_regions/data/D1p2/INALT10_FOCI/';
expe = {'FOCI1.14-II010','FOCI1.14-II011','FOCI1.14-JH027','FOCI1.19-JH037','FOCI1.19-JH039','FOCI1.14-SW128'};
addpath functions
addpath(genpath('../other_functions'));


regions = [1 2 3 4 5 6 7 8 9 10 11 12]; % repeat for area 9...
for ii= 1:length(expe)
    addpath([in_dir expe{ii}])
for aa= 1:length(regions)
    region = regions(aa); % 2 % 7 % 8
    p2_plot_TS_worker(in_dir,expe{ii},region)
end
end
