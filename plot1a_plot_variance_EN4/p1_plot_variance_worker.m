function p1_plot_variance_worker(region)

year = 1980:2019;
season = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
season_names = ['JFM';'AMJ';'JAS';'OND'];
plotdepths = [0 200 500 1000 2000];
plot_num = 1; % The plots are produced in a large grid so just going to increment plot number

%% Load settings specific to this region.
[regional_settings] = regional_settings_variance_plot(region);

%% Read global nc vars
filename = ['EN4_region' num2str(region) '.mat'];
reg = load(filename);


% Preallocate output variables
XY = nan(length(reg.lon),length(reg.lat));

S_var = nan*XY;
T_var = S_var;



%% Main loop.  Loop thru depths 

for dd = 1:length(plotdepths)
    T = [];
    S = [];
    time = [];
    depth_delete = abs(reg.depth - plotdepths(dd)); depminind = find(depth_delete == min(depth_delete)); depind = depminind(1); % X
    
    T = reg.CT(:,:,depind,:);
    S = reg.SA(:,:,depind,:);
    time = reg.time;
    
    % Generate months
    month = nan*time;
    for aa = 1:length(time)
        month(aa,1) = str2num(datestr(time(aa),'mm'));
    end

    %% Calc and plot variance
    
    % split into seasons
    for aa = 1:4
        ind = find(month == season(aa,1) | month == season(aa,2) | month == season(aa,3));
        
        T_temp = squeeze(T(:,:,:,ind));
        S_temp = squeeze(S(:,:,:,ind));
        
        % compute variance at each location
        for bb = 1:length(reg.lon)
            for cc=  1:length(reg.lat)
                S_var(bb,cc) = var(squeeze(S_temp(bb,cc,:)));
                T_var(bb,cc) = var(squeeze(T_temp(bb,cc,:)));
                
            end
        end
        
%     regional_settings.sal_contours = [10e-6 10e-5 10e-4 10e-3 10e-2];
%     regional_settings.temp_contours = [10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
%     regional_settings.Scolours = (pmkmp(4,'CubicL'));
%     regional_settings.Tcolours = (pmkmp(5,'CubicL'));
        
        %% Plot S variance
        figure(1)
        subplot(5,4,plot_num);
        set(gca,'fontsize',18);
        hold on
        pcolor(reg.lon,reg.lat,S_var'); shading flat
        caxis(regional_settings.sal_caxis);
        title(['Depth: ' (num2str(floor(reg.depth(depind)))) ' m, ' season_names(aa,:)],'fontsize',18);
        cb = colorbar;
        colormap(regional_settings.Scolours);
        set(cb,'YTick',[regional_settings.sal_labels])
        set(gca,'ColorScale','log')
        ylabel(cb,'Sal variance')
        xlim(regional_settings.xlim);
        ylim(regional_settings.ylim);
        scale = 1/cos((pi/180)*((regional_settings.ylim(1)+regional_settings.ylim(2))/2));
        daspect([scale 1 1]);
        
        %% Plot T variance
        figure(2)
        subplot(5,4,plot_num);
        set(gca,'fontsize',18);
        hold on
        % contourf(lon2d,lat2d,T_var,regional_settings.temp_contours);
        pcolor(reg.lon,reg.lat,T_var'); shading flat
        title(['Depth: ' (num2str(floor(reg.depth(depind)))) ' m, ' season_names(aa,:)],'fontsize',18);
        cb = colorbar;
        colormap(regional_settings.Tcolours);
        set(gca,'ColorScale','log')
        caxis(regional_settings.temp_caxis);
        set(cb,'YTick',[regional_settings.temp_labels])
        ylabel(cb,'Temp variance')
        xlim(regional_settings.xlim);
        ylim(regional_settings.ylim);
        scale = 1/cos((pi/180)*((regional_settings.ylim(1)+regional_settings.ylim(2))/2));
        daspect([scale 1 1]);

        plot_num = plot_num + 1;
        
    end % for each season loop
    
    disp(['Finished depth: ' num2str(reg.depth(depind)) 'm']);
end % end depth loop






%% Print completed figures

% S variance
figure(1)
width  = 2500;  % frame width
height = 2500;  % frame height
pngname = (['plots/S_variance_Region_' num2str(region) '_' regional_settings.region_name]);

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


% T variance
figure(2)
width  = 2500;  % frame width
height = 2500;  % frame height
pngname = (['plots/T_Variance_Region_' num2str(region) '_' regional_settings.region_name]);

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

close all;



