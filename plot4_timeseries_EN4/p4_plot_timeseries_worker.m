function p4_plot_timeseries_worker(region)
disp(['Working on region ' num2str(region)]);
tic
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

clear reg.sigma0

% %% NEED TO CALCULATE NEUTRAL DENSITY, now pre-calculated in EN4 files
% % Need 4D p, can't remember how to repmat higher dims
% reg.p4D = nan*reg.SA;
% for aa = 1:length(reg.time)
%     reg.p4D(:,:,:,aa) = reg.p;
% end
% % Similarly for lon,lat
% [reg.lon4D,reg.lat4D,dum,dum] = ndgrid(reg.lon,reg.lat,reg.depth,reg.time);
% 
% % Need to do some pruning otherwise 'out of memory' issues
% clear dum bathy_mask
% reg = rmfield(reg,'sigma0');
% %psal =  gsw_SP_from_SA(SA,p,long,lat)
% % psal_2D = gsw_SP_from_SA(REG.SA,0,nanmean(reg.lon),nanmean(reg.lat)); % We have salinity already in EN4
% % t = gsw_t_from_CT(SA,CT,p)
% reg.temp = gsw_t_from_CT(reg.SA,reg.CT,reg.p4D);
% % gamma_n = eos80_legacy_gamma_n(SP,t,p,long,lat)
% reg.gamma_n = eos80_legacy_gamma_n(reg.S,reg.temp,0*reg.p4D,reg.lon4D,reg.lat4D);


%% Main loop below.

if isfield(regional_settings,'Yn_thresholds') % If we have one or more density thresholds
    
    num_timeseries = length(regional_settings.Yn_thresholds(:,1));
    
    CT_timeseries = nan*ones(num_timeseries,length(reg.time));
    SA_timeseries = CT_timeseries;
    
    for ww = 1:num_timeseries % for each water mass timeseries
        
        for tt = 1:length(reg.time) % for each timestep
            % sig_temp = reg.sigma0(:,:,:,tt);
            CT_temp = reg.CT(:,:,:,tt);
            SA_temp = reg.SA(:,:,:,tt);
            Yn_temp = reg.gamma_n(:,:,:,tt);
            
            ind = find(Yn_temp >= regional_settings.Yn_thresholds(ww,1) & Yn_temp < regional_settings.Yn_thresholds(ww,2));
            
            CT_timeseries(ww,tt) = nanmean(CT_temp(ind));
            SA_timeseries(ww,tt) = nanmean(SA_temp(ind));
            
        end % end time loop
    end
    
    %elseif T or S thresholds, etc
    
end

%% Need to do some basic despiking for EN4 data...
reject_sds = 4;
for ww = 1:num_timeseries % for each water mass timeseries
    
    t_mean = nanmean(CT_timeseries(ww,:));
    s_mean = nanmean(SA_timeseries(ww,:));
    t_sd = nanstd(CT_timeseries(ww,:));
    s_sd = nanstd(SA_timeseries(ww,:));
    
    t_upper = t_mean + (t_sd * reject_sds);
    t_lower = t_mean - (t_sd * reject_sds);
    s_upper = s_mean + (s_sd * reject_sds);
    s_lower = s_mean - (s_sd * reject_sds);
    
    % remove spikes
    ind = find(CT_timeseries(ww,:) > t_upper | CT_timeseries(ww,:) < t_lower);
    CT_timeseries(ww,ind) = nan;
    SA_timeseries(ww,ind) = nan;
    
    ind = find(SA_timeseries(ww,:) > s_upper | SA_timeseries(ww,:) < s_lower);
    SA_timeseries(ww,ind) = nan;
    CT_timeseries(ww,ind) = nan;
    
    % Interpolate to fill gaps
    goodind = find(~isnan(CT_timeseries(ww,:)));
    % Vq = interp1(X,V,Xq)
    CT_timeseries(ww,:) = interp1(reg.time(goodind),CT_timeseries(ww,goodind),reg.time);
    SA_timeseries(ww,:) = interp1(reg.time(goodind),SA_timeseries(ww,goodind),reg.time);
end



%% Compute monthly anomalies
CT_seasonal = nan*ones(num_timeseries,12);
SA_seasonal = CT_seasonal;
CT_seasonal_anom = CT_seasonal;
SA_seasonal_anom = SA_seasonal;


for mm = 1:12
   ind = find(reg.month == mm);
   
   for ww = 1:num_timeseries % for each water mass timeseries
    CT_seasonal(ww,mm) = nanmean(CT_timeseries(ww,ind));
    SA_seasonal(ww,mm) = nanmean(SA_timeseries(ww,ind));
    

   end
   
end

for ww = 1:num_timeseries
    CT_seasonal_anom(ww,:) = CT_seasonal(ww,:) - mean(CT_seasonal(ww,:));
    SA_seasonal_anom(ww,:) = SA_seasonal(ww,:) - mean(SA_seasonal(ww,:));   
end

% Subtract monthly anomalies from timeseries
CT_timeseries_deseasoned = CT_timeseries*nan;
SA_timeseries_deseasoned = CT_timeseries_deseasoned;

for ww = 1:num_timeseries
    for mm = 1:12    
        ind = find(reg.month == mm);   
        CT_timeseries_deseasoned(ww,ind) = CT_timeseries(ww,ind) - CT_seasonal_anom(ww,mm);
        SA_timeseries_deseasoned(ww,ind) = SA_timeseries(ww,ind) - SA_seasonal_anom(ww,mm);
    end
end


% %% 12 month running mean
% smoothspan = 12;
% for ww = 1:num_timeseries % for each water mass timeseries
%     CT_timeseries(ww,:) = smooth(CT_timeseries(ww,:),smoothspan);
%     SA_timeseries(ww,:) = smooth(SA_timeseries(ww,:),smoothspan);
%     CT_timeseries(ww,1:smoothspan) = nan;CT_timeseries(ww,end-smoothspan:end) = nan;
%     SA_timeseries(ww,1:smoothspan) = nan;SA_timeseries(ww,end-smoothspan:end) = nan;
% end 

%% Annual means
CT_annual = nan*ones(num_timeseries,length(year));
SA_annual = CT_annual;

for ww = 1:num_timeseries % for each water mass timeseries
    for aa = 1:length(year)
        
        ind = find(reg.year == year(aa));
        CT_annual(ww,aa) = nanmean(CT_timeseries_deseasoned(ww,ind));
        SA_annual(ww,aa) = nanmean(SA_timeseries_deseasoned(ww,ind));
        
    end
end


%% Calculate trend and confidence intervals

for ww = 1:num_timeseries % for each water mass timeseries
    
    %% %%%%%%%%%%%%%% Linear regression: CT  %%%%%%%%%%%%%%%%%%%
    X = year'; tempones = ones(length(X),1); X = [X tempones]; % stats output requires a column of ones to be appended to X
    y = CT_annual(ww,:)';
    %[b,bint,r,rint,stats] = regress(y,X)    X is predictors, eg. year; y is responses, eg. salinity
    [b,bint,r,rint,stats] = regress(y,X);
    CT_lin_fit(ww,:) = b(2)+b(1)*year;
    CT_grad = b(1); CT_intercept = b(2);
    
    % Note that an anticorrelation also returns positive values of R-squared.
    %rsquared = stats(1);
    %pval = stats(3);
    
    % Get plot confidence intervals
    beta = [b(2) b(1)];
    [CT_upp_bound(ww,:),CT_low_bound(ww,:),xvals] = regression_line_ci(0.05,beta,X(:,1),y,(length(year)-1)); % for 95% confindence interval
    clear b bint r rint stats X y
    
    %% Regression analysis assuming 1D variable (dim = 1xN, Kristin code 161221)
    yfit1 = CT_lin_fit(ww,:); %this is the trend time series
    y = CT_annual(ww,:); %this is the data time series
    %t-value – find doef
    y_xc = xcorr(y-nanmean(y(:)),'coeff'); %autocorrelation
    y_xc = y_xc(ceil(length(y_xc)/2)+1:end); %only positive values
    deof_y = round(length(y)/find(y_xc<=1/exp(1),1,'first')); %deof are (number of obs)/(decorrelation or e-folding timescale)
    alphaup =  1-0.05/2; %1-0.05/2;
    t_crit_y = tinv(alphaup,deof_y); %find t_crit through matlab function
    
    yresid1 = y - yfit1;
    S2_e = (1/(deof_y-2)) * sum(yresid1.^2); % (eq. 1.33 on page 48 of internet doc)
    S_a = sqrt(S2_e / (sum((y-nanmean(y)).^2))); % (eq. 1.34 on page 48)
    % SJ NOTE: THE ABOVE COUPLE OF LINES SEEMS TO RESULT IN VALUES TOO
    % LARGE TO RESULT IN TEST PASSES. INVESTIGATE! (see https://o2.eas.gatech.edu/courses/EAS2655/week5.pdf)
    
    % significance test, t-test, 95% interval, H_0: R=0.0
    t1 = beta(2)/S_a; % (eq. 135 on page 48)
    
    % if 1: test is significant; 0: test is not significant
    if t_crit_y<abs(t1)
        CT_significant = 1;
    else
        CT_significant=0;
    end
    
    
    
    
    
    %% %%%%%%%%%%%%%% Linear regression: SA  %%%%%%%%%%%%%%%%%%%
    X = year'; tempones = ones(length(X),1); X = [X tempones]; % stats output requires a column of ones to be appended to X
    y = SA_annual(ww,:)';
    %[b,bint,r,rint,stats] = regress(y,X)    X is predictors, eg. year; y is responses, eg. salinity
    [b,bint,r,rint,stats] = regress(y,X);
    SA_lin_fit(ww,:) = b(2)+b(1)*year;
    SA_grad = b(1); SA_intercept = b(2);
    
    % Note that an anticorrelation also returns positive values of R-squared.
    %     rsquared = stats(1);
    %     pval = stats(3);
    
    % Get confidence intervals!
    beta = [b(2) b(1)];
    [SA_upp_bound(ww,:),SA_low_bound(ww,:),xvals] = regression_line_ci(0.05,beta,X(:,1),y,(length(year)-1)); % for 95% confindence interval
    clear b bint r rint stats X y
    
    %% Regression analysis assuming 1D variable (dim = 1xN, Kristin code 161221)
    yfit1 = SA_lin_fit(ww,:); %this is the trend time series
    y = SA_annual(ww,:); %this is the data time series
    %t-value – find doef
    y_xc = xcorr(y-nanmean(y(:)),'coeff'); %autocorrelation
    y_xc = y_xc(ceil(length(y_xc)/2)+1:end); %only positive values
    deof_y = round(length(y)/find(y_xc<=1/exp(1),1,'first')); %deof are (number of obs)/(decorrelation or e-folding timescale)
    alphaup =  1-0.05/2; %1-0.05/2;
    t_crit_y = tinv(alphaup,deof_y); %find t_crit through matlab function
    
    yresid1 = y - yfit1;
    S2_e = (1/(deof_y-2)) * sum(yresid1.^2); % (eq. 1.33 on page 48 of internet doc)
    S_a = sqrt(S2_e / (sum((y-nanmean(y)).^2))); % (eq. 1.34 on page 48)
    % SJ NOTE: THE ABOVE COUPLE OF LINES SEEMS TO RESULT IN VALUES TOO
    % LARGE TO RESULT IN TEST PASSES. INVESTIGATE! (see https://o2.eas.gatech.edu/courses/EAS2655/week5.pdf)
    
    % significance test, t-test, 95% interval, H_0: R=0.0
    t1 = beta(2)/S_a; % (eq. 135 on page 48)
    
    % if 1: test is significant; 0: test is not significant
    if t_crit_y<abs(t1)
        SA_significant = 1;
    else
        SA_significant=0;
    end
    
    
    
    %% Report stats to csv
    % [water mass number CT_gradient(deg/yr) CT_intercept CT_significant SA_gradient(gkg-1/yr) SA_intercept SA_significant]
    M(ww,:) = [ww CT_grad CT_intercept CT_significant SA_grad SA_intercept SA_significant];
    % M(ww,:) = [ww CT_grad CT_intercept SA_grad SA_intercept];
    
end









%% Figures and tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%% Save stats table
% dlmwrite('FILENAME',M)
dlmwrite(['plots/EN4_timeseries_REG_' num2str(region) '_' regional_settings.region_name '_STATS.csv'],M)





%% Timeseries figure
figure(1)
clf;
plot_no = 1;
% xyears = [datenum('01-01-1980') datenum('01-01-1990') datenum('01-01-2000') datenum('01-01-2010') datenum('01-01-2020')];
xyears = [1980 1990 2000 2010 2020];
for nn = 1:num_timeseries

    
    %% Tempreature timeseries
    subplot(num_timeseries,2,plot_no); 
    hold on
    
    % plot linear fit and CIs
    patchX = [year year(end:-1:1)];
    patchY = [CT_upp_bound(nn,:) CT_low_bound(nn,end:-1:1)];
    patch(patchX,patchY,[0.5 0.5 0.5],'LineStyle','none'); alpha(0.3);
    plot(year,CT_lin_fit(nn,:),'color',[0.5 0.5 0.5],'linewidth',2);
    
    %plot(reg.time,CT_timeseries(nn,:),'k');
    plot(year,CT_annual(nn,:),'k','linewidth',2);
    
    grid on;
    xticks(xyears);
    % xlim([datenum('01-01-1980') datenum('01-01-2020')]);
    xlim([1980 2020]);
    % datetick('keepticks','keeplimits');
    set(gca,'fontsize',14);
    ylabel('Cons. temp (^o C)','fontsize',14);
    title([regional_settings.WM_names{nn}])
    
    plot_no = plot_no+1;
    
    %% Salinity timeseries
    subplot(num_timeseries,2,plot_no); 
    hold on
    
    % plot linear fit and CIs
    patchX = [year year(end:-1:1)];
    patchY = [SA_upp_bound(nn,:) SA_low_bound(nn,end:-1:1)];
    patch(patchX,patchY,[0.5 0.5 0.5],'LineStyle','none'); alpha(0.3);
    plot(year,SA_lin_fit(nn,:),'color',[0.5 0.5 0.5],'linewidth',2);
    
    % plot(reg.time,SA_timeseries(nn,:),'k');
    plot(year,SA_annual(nn,:),'k','linewidth',2);
    grid on;
    xticks(xyears);
    % xlim([datenum('01-01-1980') datenum('01-01-2020')]);
    xlim([1980 2020]);
    % datetick('keepticks','keeplimits');
    set(gca,'fontsize',14);
    ylabel('Abs. Sal. (g / kg)','fontsize',14);
    title([regional_settings.WM_names{nn}])
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




% %% Seasonality figure
% figure(2)
% clf;
% plot_no = 1;
% for nn = 1:num_timeseries
% 
%     
%     %% Tempreature seasonality
%     subplot(num_timeseries,2,plot_no); 
%     hold on
%     % plot individual years in grey
%     for yy = min(reg.year):max(reg.year)
%         ind = find(reg.year == yy);   
%         plot(reg.month(ind),CT_timeseries(nn,ind),'color',[0.7 0.7 0.7]);
%     end
%     % plot monthly mean
%     plot(1:12,CT_seasonal(nn,:),'k','linewidth',2);
%     grid on;
%     xlim([1 12]);
%     set(gca,'fontsize',14);
%     xlabel('Month','fontsize',14);
%     ylabel('Cons. temp (^o C)','fontsize',14);
%     title([regional_settings.WM_names{nn}])
%     
%     plot_no = plot_no+1;
%     
%     %% Salinity seasonality
%     subplot(num_timeseries,2,plot_no); 
%     hold on
%     % plot individual years in grey
%     for yy = min(reg.year):max(reg.year)
%         ind = find(reg.year == yy);   
%         plot(reg.month(ind),SA_timeseries(nn,ind),'color',[0.7 0.7 0.7]);
%     end
%     % plot monthly mean
%     plot(1:12,SA_seasonal(nn,:),'k','linewidth',2);
%     grid on;
%     xlim([1 12]);
%     set(gca,'fontsize',14);
%     xlabel('Month','fontsize',14);
%     ylabel('Abs. Sal. (g / kg)','fontsize',14);
%     title([regional_settings.WM_names{nn}])
%     plot_no = plot_no+1;
% end
% 
% 
% 
% 
% 
% %% print figure
% width  = 2000;  % frame width
% height = 2000;  % frame height
% pngname = (['plots/seasonality/p4_Seasonality_REG_' num2str(region) '_' regional_settings.region_name]);
% 
% % set background color (outside axes)
% set(gcf,'color',[1 1 1]);
% 
% % don't change background color when printing
% set(gcf,'inverthardcopy','off');
% 
% % set size of frame to be written
% resolution=150;
% set(gcf,'paperunits','inches');
% set(gcf,'paperposition',[0 0 width height]./resolution);
% 
% % write .png file
% % the 'zbuffer' method is likely to look similar to the plot window
% print('-dpng', ['-r' num2str(resolution)], '-opengl', pngname);





toc



