function p2_plot_profiles_worker(region)

year = 1980:2019;
%season = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
%season_names = ['JFM';'AMJ';'JAS';'OND'];
basetime = datenum('01-01-1900');
% plotdepths = [0 200 500 1000 2000];
colours = (pmkmp(40,'CubicL'));

%% Read global nc vars

filename = ['1_VIKING20X.L46-KFS003_1m_19800101_19801231_grid_T_reg' num2str(region) '.nc'];
z = ncread(filename,'depth');
z2d = ncread(filename,'gdept_0'); % partial grid cells at bottom due to bathymetry
x = ncread(filename,'x');
y = ncread(filename,'y');
lon2d = ncread(filename,'glamt');
lat2d = ncread(filename,'gphit');
dx2d = ncread(filename,'e1t');
dy2d = ncread(filename,'e2t');
bathyind = ncread(filename,'mbathy');
depth = ncread(filename,'depth');
time = ncread(filename,'time');

% create real bathy values out of index bathy
bathy = nan * bathyind;
for aa = 1:length(x)
    for bb = 1:length(y)
        if bathyind(aa,bb) == 0         
        else
        bathy(aa,bb) = depth(double(bathyind(aa,bb))); 
        end
    end
end


%% Load settings specific to this region.
[regional_settings] = regional_settings_profiles_plot(region,bathy);



%% Main loop below.

% preallocate
T_mean = nan*ones(length(depth),length(year));
S_mean = T_mean;

% Loop thru years.
disp(['Creating annual means of T and S for region ' num2str(region)]);
for yy = 1:length(year) % for each year netcdf
    filename = ['1_VIKING20X.L46-KFS003_1m_' num2str(year(yy)) '0101_' num2str(year(yy)) '1231_grid_T_reg' num2str(region) '.nc'];
    % Read local nc vars
    T = ncread(filename,'votemper');
    S = ncread(filename,'vosaline');
    time = ncread(filename,'time');
    time = time + basetime;
    
    % Average ovewr time
    T = mean(T,4);
    S = mean(S,4);
    
    % Set values inside bathy from 0 to nan.
    T(T==0) = nan;
    S(S==0) = nan;
    
    %% Subset by polygon?
    if ~isempty(regional_settings.boundary_polygon)
        sz = length(lon2d(:,1)) * length(lon2d(1,:));
        lon1d = reshape(lon2d,1,sz);
        lat1d = reshape(lat2d,1,sz);
        % IN = inpolygon(X,Y,XV,YV)
        IN = inpolygon(lon1d,lat1d,regional_settings.boundary_polygon(1,:),regional_settings.boundary_polygon(2,:));
        mask = reshape(IN,length(lon2d(:,1)),length(lon2d(1,:)))';
        mask = mask';
        % Now delete everything outside mask.  Currently in loop; tidy and redo
        % with matrices...
        for aa = 1:length(depth)
                tmp = S(:,:,aa); tmp(mask==0) = nan; S(:,:,aa) = tmp;
                tmp = T(:,:,aa); tmp(mask==0) = nan; T(:,:,aa) = tmp;
        end % end depth loop
        
    end
    
    
    
    
    %% Mask using bathymetry if necessary. 
    % Set data outside bathy mask to NaN.  Can probably do this without a
    % loop but the indexing currently escapes me...
    
    for dd = 1:length(depth) % for each depth
        T_temp = T(:,:,dd);
        T_temp(regional_settings.bathy_mask==0) = nan;
        % T(:,:,dd) = T_temp;
        T_mean(dd,yy) = nanmean(nanmean(T_temp));
        
        S_temp = S(:,:,dd);
        S_temp(regional_settings.bathy_mask==0) = nan;
        % S(:,:,dd) = S_temp;
        S_mean(dd,yy) = nanmean(nanmean(S_temp));
        
    end % End of 'for each depth'
    

    
    
    
    
    
    
    disp(['Completed ' num2str(year(yy))]);
end % End year loop
    
%% Convert GSW
% Convert to TEOS-10
% [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)
SA = gsw_SA_from_SP(S_mean,depth,nanmean(nanmean(lon2d)),nanmean(nanmean(lat2d))); 
% CT = gsw_CT_from_t(SA,t,p)
CT = gsw_CT_from_t(SA,T_mean,depth);
% sigma0 = gsw_sigma0(SA,CT)
% sigma0 = gsw_sigma0(ocean.SA,ocean.CT);



    
%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



    



% Plot
figure(1)
clf

%%  temperature subplot
subplot(1,2,1);
hold on;

% plot data
for yy = 1:length(year) 
plot(CT(:,yy),depth,'linewidth',2,'color',colours(yy,:));
end


%xlabel('Absolute salinity (g/kg)','fontsize',14)
xlabel('Conservative temperature (^oC)','fontsize',14)
ylabel('Depth (m)')
set(gca,'fontsize',14);
set(gca,'ydir','reverse');
grid on;
% xlim([mins maxs]); ylim([mint maxt]);

%%  salinity subplot
subplot(1,2,2);
hold on;

% plot data
for yy = 1:length(year) 
plot(SA(:,yy),depth,'linewidth',2,'color',colours(yy,:));
end


xlabel('Absolute salinity (g/kg)','fontsize',14)
%xlabel('Conservative temperature (^oC)','fontsize',14)
ylabel('Depth (m)')
set(gca,'fontsize',14);
set(gca,'ydir','reverse');
grid on;
% xlim([mins maxs]); ylim([mint maxt]);





% Colorbar
colormap(colours);
caxis([1979.5 2019.5]);
cb = colorbar;
set(cb,'fontsize',14);


%% print figure
width  = 2000;  % frame width
height = 1500;  % frame height
pngname = (['plots/p2_V20_profiles_coloured_time_REG_' num2str(region) '_' regional_settings.region_name]);

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
