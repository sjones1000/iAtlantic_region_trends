% EN4 regions map
clear; close all;

addpath(genpath('G:\Work_computer_sync\MATLAB_functions'));
addpath(genpath('H:\EN4'));


minlon = -80;
maxlon = 22;
minlat = -42;
maxlat = 75;




%% Plot map
figure(1); hold on;
m_proj('mercator','longitudes',[minlon maxlon],'latitudes',[minlat maxlat]);

% bathymetry
bath = load('GEBCO_world_1D.mat');
bath.bathy = bath.bathy';
ind = find(bath.longitude >= minlon & bath.longitude < maxlon); bath.longitude = bath.longitude(ind); bath.bathy = bath.bathy(:,ind);
ind = find(bath.latitude >= minlat & bath.latitude < maxlat); bath.latitude = bath.latitude(ind); bath.bathy = bath.bathy(ind,:);
bathinv = bath.bathy.*-1;
m_contour(bath.longitude,bath.latitude,bathinv,[0 0],'color','k','linewidth',2);
m_contour(bath.longitude,bath.latitude,bathinv,[2000 4000 6000],'color',[0.7 0.7 0.7],'linewidth',1.5);

%% Load and plot surface data from each EN4 region

for aa = 1:12
filename = ['EN4_region' num2str(aa) '.mat'];
reg = load(filename);

meanT = mean(reg.CT(:,:,1,:),4);

textx = mean(reg.lon);
texty = mean(reg.lat);
m_text(textx,texty,num2str(aa),'fontsize',20);

m_pcolor(reg.lon,reg.lat,meanT'); shading flat

clear reg meanT

end

caxis([0 25]);
colormap(pmkmp(25,'CubicL'));
set(gca,'fontsize',20);
cb = colorbar;
ylabel(cb,'Annual mean surface temperature (^o C)','fontsize',20);


% tidying uup
m_grid('box','fancy','fontsize',20,'color','k','tickdir','in');




%% print figure
width  = 2000;  % frame width
height = 2000;  % frame height
pngname = ('plots/p3_EN4_region_map');

% set background color (outside axes)
set(gcf,'color',[1 1 1]);

% don't change background color when printing
set(gcf,'inverthardcopy','off');

% set size of frame to be written
resolution=150;
set(gcf,'paperunits','inches');
set(gcf,'paperposition',[0 0 width height]./resolution);

% write .png file
% the 'zbuffer' method is likely to look similar to the plot window
print('-dpng', ['-r' num2str(resolution)], '-opengl', pngname);