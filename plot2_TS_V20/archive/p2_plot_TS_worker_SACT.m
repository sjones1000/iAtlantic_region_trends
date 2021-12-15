function p2_plot_TS_worker(region)

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
[regional_settings] = regional_settings_TS_plot(region,bathy);



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

%% Generate density contours
mint = min(min(CT))-0.2; maxt = max(max(CT))+0.2;
mins = min(min(SA))-0.05; maxs = max(max(SA))+0.05;
xdim = round((maxs-mins)./0.005+1);
ydim=round((maxt-mint)./0.05+1);
sig0=zeros(ydim,xdim);
gamma_n=zeros(ydim,xdim);
sig0_y=((1:ydim)-1)*0.05+mint;
sig0_x=((1:xdim)-1)*0.005+mins;

% Conversion to neutral density for contours
% sig0_x is salinity increments, y is temp increments

[SA_2D,CT_2D] = meshgrid(sig0_x,sig0_y);
%repmat(nanmean(sig0_x),1,length(sig0_y));

%psal =  gsw_SP_from_SA(SA,p,long,lat)
psal_2D = gsw_SP_from_SA(SA_2D,0,nanmean(nanmean(lon2d)),nanmean(nanmean(lat2d)));
% t = gsw_t_from_CT(SA,CT,p)
temp_2D = gsw_t_from_CT(SA_2D,CT_2D,0);
% gamma_n = eos80_legacy_gamma_n(SP,t,p,long,lat)
gamma_n = eos80_legacy_gamma_n(psal_2D,temp_2D,0,nanmean(nanmean(lon2d)),nanmean(nanmean(lat2d)));
% sigma0 = gsw_sigma0(SA,CT)
sig0 = gsw_sigma0(SA_2D,CT_2D);


    



% Plot
figure(1)
clf
hold on;
mindens = floor(min(min(sig0)));
maxdens = ceil(max(max(sig0)));
% levels = mindens:0.25:maxdens;
% [c,h]=contour(sig0_x,sig0_y,sig0,levels,'color',[0.4 0.4 0.4]);
% clabel(c,h,'LabelSpacing',1000,'color',[0.4 0.4 0.4],'fontsize',14);

Ynlevels = mindens:0.1:maxdens;
[d,j]=contour(sig0_x,sig0_y,gamma_n,Ynlevels,'color',[0.4 0.4 0.7]);
clabel(d,j,'LabelSpacing',1000,'color',[0.4 0.4 0.4],'fontsize',14);


xlabel('Absolute salinity (g/kg)','fontsize',20)
ylabel('Conservative temperature (^oC)','fontsize',20)
set(gca,'fontsize',20);
% xlim([mins maxs]); ylim([mint maxt]);

% plot data
for yy = 1:length(year) 
plot(SA(:,yy),CT(:,yy),'linewidth',2,'color',colours(yy,:));
end

% plot water masses
WM = load('source_WM.mat');

for aa = regional_settings.water_mass_plot
    if ~isnan(WM.CT(aa,2)) % if there is an upper and lower line
        line([WM.SA(aa,1) WM.SA(aa,2)],[WM.CT(aa,1) WM.CT(aa,2)],'color','k','linewidth',2);
        plot(WM.SA(aa,1),WM.CT(aa,1),'ok','linewidth',2);
        plot(WM.SA(aa,2),WM.CT(aa,2),'ok','linewidth',2);
        text(((WM.SA(aa,1)+WM.SA(aa,2))/2),((WM.CT(aa,1)+WM.CT(aa,2))/2),WM.Name(aa),'fontsize',14);
    else % if there is just a point
        plot(WM.SA(aa,1),WM.CT(aa,1),'+k','linewidth',2);
        text(WM.SA(aa,1),WM.CT(aa,1),WM.Name(aa),'fontsize',14);
    end
end

% Colorbar
colormap(colours);
caxis([1979.5 2019.5]);
cb = colorbar;
set(cb,'fontsize',20);

% Axis limits
if ~isempty(regional_settings.axis_limits)
    % [x1 x2 y1 y2];
    xlim([regional_settings.axis_limits(1) regional_settings.axis_limits(2)]);
    ylim([regional_settings.axis_limits(3) regional_settings.axis_limits(4)]);
else
    xlim([mins maxs]);
    ylim([mint maxt]);
end

%% print figure
width  = 1500;  % frame width
height = 1500;  % frame height
pngname = (['plots/p2_V20_TS_coloured_time_REG_' num2str(region) '_' regional_settings.region_name]);

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
