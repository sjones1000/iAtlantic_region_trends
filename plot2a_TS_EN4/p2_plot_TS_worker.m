function p2_plot_TS_worker(region)

year = 1980:2020;
%season = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
%season_names = ['JFM';'AMJ';'JAS';'OND'];
% plotdepths = [0 200 500 1000 2000];
colours = (pmkmp(41,'CubicL'));

%% Read global nc vars
% filename = ['1_VIKING20X.L46-KFS003_1m_19800101_19801231_grid_T_reg' num2str(region) '.nc'];
% z = ncread(filename,'depth');
% z2d = ncread(filename,'gdept_0'); % partial grid cells at bottom due to bathymetry
% x = ncread(filename,'x');
% y = ncread(filename,'y');
% lon2d = ncread(filename,'glamt');
% lat2d = ncread(filename,'gphit');
% dx2d = ncread(filename,'e1t');
% dy2d = ncread(filename,'e2t');
% bathyind = ncread(filename,'mbathy');
% depth = ncread(filename,'depth');
% time = ncread(filename,'time');

%% Read global nc vars
filename = ['EN4_region' num2str(region) '.mat'];
reg = load(filename);

% Create year variable
reg.year = nan*reg.time;
for aa = 1:length(reg.time)
   reg.year(aa,1) = str2num(datestr(reg.time(aa),'yyyy'));
end

% back-calculate in-situ T from ptemp
p4D = nan*reg.SA;
for aa = 1:length(reg.time)
    p4D(:,:,:,aa) = reg.p;
end
% t = gsw_t_from_pt0(SA,pt0,p)
reg.T = gsw_t_from_pt0(reg.SA,reg.Tptemp,p4D);
clear p4D

%% Compute bathy associated with each EN4 cell.
bath = load('G:\Work_computer_sync\MATLAB_functions\m_map\GEBCO_world_1D.mat');
% bath.bathy = bath.bathy';
% ind = find(bath.longitude >= minlon & bath.longitude < maxlon); bath.longitude = bath.longitude(ind); bath.bathy = bath.bathy(:,ind);
% ind = find(bath.latitude >= minlat & bath.latitude < maxlat); bath.latitude = bath.latitude(ind); bath.bathy = bath.bathy(ind,:);


EN4_bathy = nan(length(reg.lon),length(reg.lat));

for aa = 1:length(reg.lon)
   for bb = 1:length(reg.lat)
       lon_delete = abs(reg.lon(aa) - bath.longitude); lonminind = find(lon_delete == min(lon_delete)); lonind = lonminind(1);
       lat_delete = abs(reg.lat(bb) - bath.latitude); latminind = find(lat_delete == min(lat_delete)); latind = latminind(1);
       EN4_bathy(aa,bb) = bath.bathy(lonind,latind);
   end
end

EN4_bathy = EN4_bathy.*-1;

% % create real bathy values out of index bathy
% bathy = nan * bathyind;
% for aa = 1:length(x)
%     for bb = 1:length(y)
%         if bathyind(aa,bb) == 0         
%         else
%         bathy(aa,bb) = depth(double(bathyind(aa,bb))); 
%         end
%     end
% end


%% Load settings specific to this region.
[regional_settings] = regional_settings_TS_plot(region,EN4_bathy);


%% Subset by polygon?
if ~isempty(regional_settings.boundary_polygon)
    sz = length(reg.lon) * length(reg.lat);
    [lonm,latm] = meshgrid(reg.lon,reg.lat);
    lon1d = reshape(lonm,1,sz);
    lat1d = reshape(latm,1,sz);
    % IN = inpolygon(X,Y,XV,YV)
    IN = inpolygon(lon1d,lat1d,regional_settings.boundary_polygon(1,:),regional_settings.boundary_polygon(2,:));
    mask = reshape(IN,length(reg.lat),length(reg.lon))';

    % Now delete everything outside mask.  Currently in loop; tidy and redo
    % with matrices...
    for aa = 1:length(reg.depth)
       for bb = 1:length(reg.time)
           tmp = reg.CT(:,:,aa,bb); tmp(mask==0) = nan; reg.CT(:,:,aa,bb) = tmp;
           tmp = reg.SA(:,:,aa,bb); tmp(mask==0) = nan; reg.SA(:,:,aa,bb) = tmp;
           tmp = reg.S(:,:,aa,bb); tmp(mask==0) = nan; reg.S(:,:,aa,bb) = tmp;
           tmp = reg.T(:,:,aa,bb); tmp(mask==0) = nan; reg.T(:,:,aa,bb) = tmp;
           tmp = reg.Tptemp(:,:,aa,bb); tmp(mask==0) = nan; reg.Tptemp(:,:,aa,bb) = tmp;
           tmp = reg.sigma0(:,:,aa,bb); tmp(mask==0) = nan; reg.sigma0(:,:,aa,bb) = tmp;
       end % end time loop
    end % end depth loop

end





%% Main loop below.

% preallocate
T_mean = nan*ones(length(reg.depth),length(year));
S_mean = T_mean;

% Loop thru years.
disp(['Creating annual means of T and S for region ' num2str(region)]);
for yy = 1:length(year) % for each year netcdf
    
    ind = find(reg.year == year(yy));
    
    % Read local nc vars
    T = reg.T(:,:,:,ind);
    S = reg.S(:,:,:,ind);
    time = reg.time(ind);
    
    % Average ovewr time
    T = mean(T,4);
    S = mean(S,4);
    
    %% Mask using bathymetry if necessary. 
    % Set data outside bathy mask to NaN.  Can probably do this without a
    % loop but the indexing currently escapes me...
    
    for dd = 1:length(reg.depth) % for each depth
        T_temp = T(:,:,dd);
        T_temp(regional_settings.bathy_mask==0) = nan;
        % T(:,:,dd) = T_temp;
        T_mean(dd,yy) = nanmean(nanmean(T_temp));
        
        S_temp = S(:,:,dd);
        S_temp(regional_settings.bathy_mask==0) = nan;
        % S(:,:,dd) = S_temp;
        S_mean(dd,yy) = nanmean(nanmean(S_temp));
        
    end % End of 'for each depth'

end % End year loop
    
% Replace 3 dimensional S and T with averaged versions
T = T_mean;
S = S_mean;



    
%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%% Generate density contours
mint = min(min(T))-0.2; maxt = max(max(T))+0.2;
mins = min(min(S))-0.05; maxs = max(max(S))+0.05;
xdim = linspace(mins,maxs,100);
ydim = linspace(mint,maxt,100);
[S_2D,T_2D] = meshgrid(xdim,ydim);
% gamma_n = eos80_legacy_gamma_n(SP,t,p,long,lat)
gamma_n = eos80_legacy_gamma_n(S_2D,T_2D,0,nanmean(reg.lon),nanmean(reg.lat));







% Plot
figure(1)
clf
hold on;
mindens = floor(min(min(gamma_n)));
maxdens = ceil(max(max(gamma_n)));
Ynlevels = mindens:0.1:maxdens;
[d,j]=contour(xdim,ydim,gamma_n,Ynlevels,'color',[0.4 0.4 0.7]);
clabel(d,j,'LabelSpacing',1000,'color',[0.4 0.4 0.4],'fontsize',14);


xlabel('Practical salinity','fontsize',20)
ylabel('Temperature (^oC)','fontsize',20)
set(gca,'fontsize',20);
% xlim([mins maxs]); ylim([mint maxt]);

% plot data
for yy = 1:length(year) 
plot(S(:,yy),T(:,yy),'linewidth',2,'color',colours(yy,:));
end

% % plot water masses
% WM = load('source_WM.mat');
% 
% for aa = regional_settings.water_mass_plot
%     if ~isnan(WM.CT(aa,2)) % if there is an upper and lower line
%         line([WM.SA(aa,1) WM.SA(aa,2)],[WM.CT(aa,1) WM.CT(aa,2)],'color','k','linewidth',2);
%         plot(WM.SA(aa,1),WM.CT(aa,1),'ok','linewidth',2);
%         plot(WM.SA(aa,2),WM.CT(aa,2),'ok','linewidth',2);
%         text(((WM.SA(aa,1)+WM.SA(aa,2))/2),((WM.CT(aa,1)+WM.CT(aa,2))/2),WM.Name(aa),'fontsize',14);
%     else % if there is just a point
%         plot(WM.SA(aa,1),WM.CT(aa,1),'+k','linewidth',2);
%         text(WM.SA(aa,1),WM.CT(aa,1),WM.Name(aa),'fontsize',14);
%     end
% end

% Colorbar
colormap(colours);
caxis([1979.5 2020.5]);
cb = colorbar;
set(cb,'fontsize',20);

% Axis limits
if region == 4
    xlim([34.5 maxs]);
    ylim([mint maxt]);
else
    xlim([mins maxs]);
    ylim([mint maxt]);
end
% if ~isempty(regional_settings.axis_limits)
%     % [x1 x2 y1 y2];
%     xlim([regional_settings.axis_limits(1) regional_settings.axis_limits(2)]);
%     ylim([regional_settings.axis_limits(3) regional_settings.axis_limits(4)]);
% else
%     xlim([mins maxs]);
%     ylim([mint maxt]);
% end

%% print figure
width  = 1500;  % frame width
height = 1500;  % frame height
pngname = (['plots/t_psal/p2_TS_coloured_time_REG_' num2str(region) '_' regional_settings.region_name]);

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

