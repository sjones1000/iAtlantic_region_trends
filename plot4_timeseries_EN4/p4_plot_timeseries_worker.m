function p4_plot_timeseries_worker(region)

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

% Create year / month variable
reg.year = nan*reg.time;
for aa = 1:length(reg.time)
   reg.year(aa,1) = str2num(datestr(reg.time(aa),'yyyy'));
   reg.month(aa,1) = str2num(datestr(reg.time(aa),'mm'));
end

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
[regional_settings] = regional_settings_timeseries_plot(region,EN4_bathy);





% preallocate
T_mean = nan*ones(length(reg.depth),length(year));
S_mean = T_mean;

% Mask bathy if required
bathy_mask = repmat(regional_settings.bathy_mask,[1,1,length(reg.depth),length(reg.time)]);
reg.CT(bathy_mask==0) = nan;
reg.SA(bathy_mask==0) = nan;
reg.sigma0(bathy_mask==0) = nan;

%% Main loop below.

if isfield(regional_settings,'sig0_thresholds') % If we have one or more density thresholds
    
    num_timeseries = length(regional_settings.sig0_thresholds(:,1));
    
    CT_timeseries = nan*ones(num_timeseries,length(reg.time));
    SA_timeseries = CT_timeseries;
    
    for ww = 1:num_timeseries % for each water mass timeseries
        
        for tt = 1:length(reg.time) % for each timestep
            sig_temp = reg.sigma0(:,:,:,tt);
            CT_temp = reg.CT(:,:,:,tt);
            SA_temp = reg.SA(:,:,:,tt);
            
            ind = find(sig_temp >= regional_settings.sig0_thresholds(ww,1) & sig_temp < regional_settings.sig0_thresholds(ww,2));
            
            CT_timeseries(ww,tt) = nanmean(CT_temp(ind));
            SA_timeseries(ww,tt) = nanmean(SA_temp(ind));
            
        end % end time loop
    end
    
    %elseif T or S thresholds, etc
    
end

%% Also break down by month

CT_seasonal = nan*ones(num_timeseries,12);
SA_seasonal = CT_seasonal;

for mm = 1:12
   ind = find(reg.month == mm);
   
   for ww = 1:num_timeseries % for each water mass timeseries
    CT_seasonal(ww,mm) = nanmean(CT_timeseries(ww,ind));
    SA_seasonal(ww,mm) = nanmean(SA_timeseries(ww,ind));
   end
   
end







    
%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%% Timeseries figure
figure(1)
clf;
plot_no = 1;
for nn = 1:num_timeseries

    
    %% Tempreature timeseries
    subplot(3,2,plot_no); 
    plot(reg.time,CT_timeseries(nn,:),'k');
    grid on;
    xlim([datenum(reg.time(1)) datenum(reg.time(end))]);
    datetick('keepticks','keeplimits');
    set(gca,'fontsize',14);
    ylabel('Cons. temp (^o C)','fontsize',14);
    title(['Temperature: ' regional_settings.WM_names(nn)])
    
    plot_no = plot_no+1;
    
    %% Salinity timeseries
    subplot(3,2,plot_no); 
    plot(reg.time,SA_timeseries(nn,:),'k');
    grid on;
    xlim([datenum(reg.time(1)) datenum(reg.time(end))]);
    datetick('keepticks','keeplimits');
    set(gca,'fontsize',14);
    ylabel('Abs. Sal. (g / kg)','fontsize',14);
    title(['Salinity: ' regional_settings.WM_names(nn)])
    plot_no = plot_no+1;
end





%% print figure
width  = 2000;  % frame width
height = 2000;  % frame height
pngname = (['plots/p4_timeseries_REG_' num2str(region) '_' regional_settings.region_name]);

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




%% Seasonality figure
figure(2)
clf;
plot_no = 1;
for nn = 1:num_timeseries

    
    %% Tempreature seasonality
    subplot(3,2,plot_no); 
    hold on
    % plot individual years in grey
    for yy = min(reg.year):max(reg.year)
        ind = find(reg.year == yy);   
        plot(reg.month(ind),CT_timeseries(nn,ind),'color',[0.7 0.7 0.7]);
    end
    % plot monthly mean
    plot(1:12,CT_seasonal(nn,:),'k','linewidth',2);
    grid on;
    xlim([1 12]);
    set(gca,'fontsize',14);
    xlabel('Month','fontsize',14);
    ylabel('Cons. temp (^o C)','fontsize',14);
    title(['Temperature: ' regional_settings.WM_names(nn)])
    
    plot_no = plot_no+1;
    
    %% Salinity seasonality
    subplot(3,2,plot_no); 
    hold on
    % plot individual years in grey
    for yy = min(reg.year):max(reg.year)
        ind = find(reg.year == yy);   
        plot(reg.month(ind),SA_timeseries(nn,ind),'color',[0.7 0.7 0.7]);
    end
    % plot monthly mean
    plot(1:12,SA_seasonal(nn,:),'k','linewidth',2);
    grid on;
    xlim([1 12]);
    set(gca,'fontsize',14);
    xlabel('Month','fontsize',14);
    ylabel('Abs. Sal. (g / kg)','fontsize',14);
    title(['Salinity: ' regional_settings.WM_names(nn)])
    plot_no = plot_no+1;
end





%% print figure
width  = 2000;  % frame width
height = 2000;  % frame height
pngname = (['plots/p4_Seasonality_REG_' num2str(region) '_' regional_settings.region_name]);

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




