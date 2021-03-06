function p4_plot_timeseries_worker(path,region)
disp(['Working on region ' num2str(region)]);
tic

expe = {'FOCI1.14-II010','FOCI1.14-II011','FOCI1.14-JH027','FOCI1.19-JH037','FOCI1.19-JH039','FOCI1.14-SW128'};
eN = length(expe);
year = 2014:2070;
basetime = datenum('01-01-1900');

exp_idx=1;
addpath([path expe{exp_idx}])
filelist= dir([path expe{exp_idx} '/' '*_reg' num2str(region) '.nc']);
files_gn=dir([path,expe{exp_idx},'/gamma_n/','*_reg' num2str(region) '.nc']);


%% Read global nc vars
filename = filelist(1).name;
z = ncread(filename,'depth');
z2d = ncread(filename,'gdept_0'); % partial grid cells at bottom due to bathymetry
x = ncread(filename,'x');
y = ncread(filename,'y');
lon2d = ncread(filename,'glamt');
lat2d = ncread(filename,'gphit');
tmask = ncread(filename,'tmask');
tmask = repmat(tmask,1,1,1,12);
dx2d = ncread(filename,'e1t');
dy2d = ncread(filename,'e2t');
bathyind = ncread(filename,'mbathy');
depth = ncread(filename,'depth');
time = ncread(filename,'time');

% Need a 3d depth var for later
depth3d = repmat(depth,1,length(lon2d(:,1)),length(lon2d(1,:)));
depth3d = permute(depth3d,[2 3 1]);
depth4d = repmat(depth3d,1,1,1,12);

% Same for lon and lat
lon3d = repmat(lon2d,1,1,length(depth));
lon4d = repmat(lon3d,1,1,1,12);
lat3d = repmat(lat2d,1,1,length(depth));
lat4d = repmat(lat3d,1,1,1,12);

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

% Create month / year, global time variable
inc = 1;
for aa = 1:length(year)
    for bb = 1:12
        reg.year(1,inc) = year(aa);
        reg.month(1,inc) = bb;
        reg.time(1,inc) = datenum([num2str(bb) '-15-' num2str(year(aa))]);
        inc = inc+1;
    end
end

%% Load settings specific to this region.
[regional_settings] = regional_settings_timeseries_plot(region,bathy);
num_timeseries = length(regional_settings.Yn_thresholds(:,1));

CT_annual_all = nan*ones(eN,num_timeseries,length(year));
SA_annual_all = CT_annual_all;


%% Main loop below. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for exp_idx=1:eN
    disp(expe{exp_idx})
    addpath([path expe{exp_idx}])
    filelist= dir([path expe{exp_idx} '/' '*_reg' num2str(region) '.nc']);
    files_gn=dir([path,expe{exp_idx},'/gamma_n/','*_reg' num2str(region) '.nc']);

    % Create the empty variables to take the full timeseries
    CT_timeseries = nan*ones(num_timeseries,(length(time)*length(year)));
    SA_timeseries = CT_timeseries;


    global_time_start = 0;
    for yy = 1:length(filelist) % for each year netcdf
        filename = filelist(yy).name;
        % Read local nc vars
        T = ncread(filename,'votemper');
        S = ncread(filename,'vosaline');
        time = ncread(filename,'time');
        time = time + basetime;

        % Set values inside bathy from 0 to nan.
        T(tmask==0) = nan;
        S(tmask==0) = nan;

        Yn_filename = [path,expe{exp_idx},'/gamma_n/',files_gn(yy).name];
        gamma_n = ncread(Yn_filename,'gamma_n');

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
                for dd = 1:length(depth)
                    for tt = 1:length(time)
                        tmp = S(:,:,dd,tt); tmp(mask==0) = nan; S(:,:,dd,tt) = tmp;
                        tmp = T(:,:,dd,tt); tmp(mask==0) = nan; T(:,:,dd,tt) = tmp;
                    end
                end % end depth loop
            end



        %% Mask using bathymetry if necessary.
        % Set data outside bathy mask to NaN.  Can probably do this without a
        % loop but the indexing currently escapes me...

        for dd = 1:length(depth) % for each depth
            for tt = 1:length(time) % for each time step
                T_temp = T(:,:,dd,tt);
                T_temp(regional_settings.bathy_mask==0) = nan;
                T(:,:,dd,tt) = T_temp;

                S_temp = S(:,:,dd,tt);
                S_temp(regional_settings.bathy_mask==0) = nan;
                S(:,:,dd,tt) = S_temp;
            end
        end % End of 'for each depth'






        % Convert GSW - Note this is just for one timestep
        % NOTE THE GSW TASKS ARE BY FAR THE SLOWEST PART OF THIS SCRIPT.
        % MAYBE THERE'S A WAY TO BYPASS UNTIL THE AVERAGE IS COMPLETED...
        % [SA, in_ocean] = gsw_SA_from_SP(SP,p,long,lat)
        % SA = gsw_SA_from_SP(S,depth4d,nanmean(nanmean(lon2d)),nanmean(nanmean(lat2d)));
        SA = gsw_SA_from_SP(S,depth4d,lon4d,lat4d);
        % CT = gsw_CT_from_t(SA,t,p)
        CT = gsw_CT_from_t(SA,T,depth4d);

        for ww = 1:num_timeseries % for each water mass timeseries

            for tt = 1:12 % for each timestep
                CT_temp = CT(:,:,:,tt);
                SA_temp = SA(:,:,:,tt);
                Yn_temp = gamma_n(:,:,:,tt);

                ind = find(Yn_temp >= regional_settings.Yn_thresholds(ww,1) & Yn_temp < regional_settings.Yn_thresholds(ww,2));

                CT_timeseries(ww,tt+global_time_start) = nanmean(CT_temp(ind));
                SA_timeseries(ww,tt+global_time_start) = nanmean(SA_temp(ind));
                % tt % print loop number
            end % end time loop
        end

        global_time_start = global_time_start + 12;
        disp(['Finished year ' num2str(year(yy))]);
    end % End of 'for each year (V20 file)' loop


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
    CT_annual_all(exp_idx,:,:) = CT_annual;
    SA_annual_all(exp_idx,:,:) = SA_annual;

end

%% Calculate trend and confidence intervals
CT_annual = squeeze(nanmean(CT_annual_all,1));
SA_annual = squeeze(nanmean(SA_annual_all,1));

for ww = 1:num_timeseries % for each water mass timeseries
    if sum(~isnan(CT_annual(ww,:)))>0
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
    %t-value ??? find doef
    y_xc = xcorr(y-nanmean(y(:)),'coeff'); %autocorrelation
    y_xc = y_xc(ceil(length(y_xc)/2)+1:end); %only positive values
    deof_y = round(length(y)/find(y_xc<=1/exp(1),1,'first')); %deof are (number of obs)/(decorrelation or e-folding timescale)
    alphaup =  1-0.05/2; %1-0.05/2;
    t_crit_y = tinv(alphaup,deof_y); %find t_crit through matlab function
    
    yresid = y - yfit1;
    yresid1 = sum(yresid.^2); % sum of the squared residuals
    ytotal1 = (length(y)-1) * var(y); % n * variance
    
    % adjusted coefficient of determination
    rsq_y = 1 - yresid1/ytotal1*(length(y)-1)/(length(y)-2);
    %rsq_y = rsquared;
    
    % significance test, t-test, 95% interval, H_0: R=0.0
    t1 = sqrt((rsq_y*(deof_y-2))/(1-rsq_y));
    if t_crit_y<abs(t1)
        CT_significant = 1;
    else
        CT_significant = 0;
    end
    
    clear y_xc deof_y alphaup t_crit_y yresid yresid1 ytotal1 rsq_y t1

    
    
    
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
    %t-value ??? find doef
    y_xc = xcorr(y-nanmean(y(:)),'coeff'); %autocorrelation
    y_xc = y_xc(ceil(length(y_xc)/2)+1:end); %only positive values
    deof_y = round(length(y)/find(y_xc<=1/exp(1),1,'first')); %deof are (number of obs)/(decorrelation or e-folding timescale)
    alphaup =  1-0.05/2; %1-0.05/2;
    t_crit_y = tinv(alphaup,deof_y); %find t_crit through matlab function
    
    yresid = y - yfit1;
    yresid1 = sum(yresid.^2); % sum of the squared residuals
    ytotal1 = (length(y)-1) * var(y); % n * variance
    
    % adjusted coefficient of determination
    rsq_y = 1 - yresid1/ytotal1*(length(y)-1)/(length(y)-2);
    %rsq_y = rsquared;    
    
    % significance test, t-test, 95% interval, H_0: R=0.0
    t1 = sqrt((rsq_y*(deof_y-2))/(1-rsq_y));
    if t_crit_y<abs(t1)
        SA_significant = 1;
    else
        SA_significant = 0;
    end
    
    clear y_xc deof_y alphaup t_crit_y yresid yresid1 ytotal1 rsq_y t1
    
    %% Report stats to csv
    % [water mass number CT_gradient(deg/yr) CT_intercept CT_significant SA_gradient(gkg-1/yr) SA_intercept SA_significant]
    M(ww,:) = [ww CT_grad.*10 CT_intercept CT_significant SA_grad.*10 SA_intercept SA_significant];
    end
end









%% Figures and tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%% Save stats table
% dlmwrite('FILENAME',M)
dlmwrite(['plots/INALT20_timeseries_REG_' num2str(region) '_' regional_settings.region_name '_STATS.csv'],M)





%% Timeseries figure
figure(1)
clf;
plot_no = 1;
% xyears = [datenum('01-01-1980') datenum('01-01-1990') datenum('01-01-2000') datenum('01-01-2010') datenum('01-01-2020')];
xyears = [2020 2030 2040 2050 2060 2070];
for nn = 1:num_timeseries

    
    %% Temperature timeseries
    subplot(num_timeseries,2,plot_no); 
    hold on
    
    %plot(reg.time,CT_timeseries(nn,:),'k');
    plot(year,squeeze(CT_annual_all(:,nn,:)),'color',[.5 .5 .5],'linewidth',1)
    plot(year,CT_annual(nn,:),'k','linewidth',2);
    
    % plot linear fit and CIs
    patchX = [year year(end:-1:1)];
    patchY = [CT_upp_bound(nn,:) CT_low_bound(nn,end:-1:1)];
    patch(patchX,patchY,[0.6350 0.0780 0.1840],'LineStyle','none'); alpha(0.7);
    plot(year,CT_lin_fit(nn,:),'color',[0.6350 0.0780 0.1840],'linewidth',2);alpha(0.5);
    grid on;
    xticks(xyears);
    
    % xlim([datenum('01-01-1980') datenum('01-01-2020')]);
    xlim([year(1) year(end)]);
    % datetick('keepticks','keeplimits');
    set(gca,'fontsize',14);
    ylabel('Cons. temp (^o C)','fontsize',14);
    title([regional_settings.WM_names{nn}])
    
    plot_no = plot_no+1;
    
    %% Salinity timeseries
    subplot(num_timeseries,2,plot_no); 
    hold on
    
    % plot(reg.time,SA_timeseries(nn,:),'k');
    plot(year,squeeze(SA_annual_all(:,nn,:)),'color',[.5 .5 .5],'linewidth',1)
    plot(year,SA_annual(nn,:),'k','linewidth',2);
    
    % plot linear fit and CIs
    patchX = [year year(end:-1:1)];
    patchY = [SA_upp_bound(nn,:) SA_low_bound(nn,end:-1:1)];
    patch(patchX,patchY,[0.6350 0.0780 0.1840],'LineStyle','none'); alpha(0.7);
    plot(year,SA_lin_fit(nn,:),'color',[0.6350 0.0780 0.1840],'linewidth',2);alpha(0.5);
    
    grid on;
    xticks(xyears);
    % xlim([datenum('01-01-1980') datenum('01-01-2020')]);
    xlim([year(1) year(end)]);
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
