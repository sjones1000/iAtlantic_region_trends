function p1_calc_Yn_worker(region)


disp(['************** Working on Region ' num2str(region) ' *******************']);
disp(['******************************************************']);
year = 1980:2019;


%% Read global nc vars
filename = ['1_VIKING20X.L46-KFS003_1m_19800101_19801231_grid_T_reg' num2str(region) '.nc'];
z = ncread(filename,'depth');
z2d = ncread(filename,'gdept_0'); % partial grid cells at bottom due to bathymetry
x = double(ncread(filename,'x'));
y = double(ncread(filename,'y'));
lon2d = ncread(filename,'glamt');
lat2d = ncread(filename,'gphit');
dx2d = ncread(filename,'e1t');
dy2d = ncread(filename,'e2t');
bathyind = ncread(filename,'mbathy');
depth = ncread(filename,'depth');
time = ncread(filename,'time');

depth = depth.*-1;

% Need a 3d depth var
depth3d = repmat(depth,1,length(lon2d(:,1)),length(lon2d(1,:)));
depth3d = permute(depth3d,[2 3 1]);
% depth4d = repmat(depth3d,1,1,1,12);

lon3d = repmat(lon2d,1,1,length(depth));
lat3d = repmat(lat2d,1,1,length(depth));

% Calc pressure
% p = gsw_p_from_z(z,lat);
p3d = gsw_p_from_z(depth3d,lat3d);



for yy = 1:length(year) % for each year netcdf
    tic
    
    % Preallocate empty gamma_n
    gamma_n = nan*ones(length(x),length(y),length(depth),length(time));
    
    % Load T and S for this year
    filename = ['1_VIKING20X.L46-KFS003_1m_' num2str(year(yy)) '0101_' num2str(year(yy)) '1231_grid_T_reg' num2str(region) '.nc'];
    T = ncread(filename,'votemper');
    S = ncread(filename,'vosaline');
    
    % Set values inside bathy from 0 to nan.
    T(T==0) = nan;
    S(S==0) = nan;
    
    
    
    
    
    % Work on one month at a time
    for mm = 1:12
        % Calc annual mean
        T3d = squeeze(T(:,:,:,mm));
        S3d = squeeze(S(:,:,:,mm));
        
        %% Calc Yn!
        % gamma_n = eos80_legacy_gamma_n(SP,t,p,long,lat)
        gamma_n(:,:,:,mm) = eos80_legacy_gamma_n(S3d,T3d,p3d,lon3d,lat3d);
        
        disp(['Month ' num2str(mm) ' *Comment out when not needed*']); 
    end % end of 'for each month' loop
    
    
    
    
    
    %% Write netcdf
    fileout = ['I:/iAtlantic_overflow/VIKING20/D1p2/gamma_n/V20X_region_' num2str(region) '_' num2str(year(yy)) '.nc']; % <- Save data location
    ncid = netcdf.create(fileout,'NETCDF4'); % open ncfile
    
    % Create dimensions
    dimid_x = netcdf.defDim(ncid,'x',length(x));
    dimid_y = netcdf.defDim(ncid,'y',length(y));
    dimid_z = netcdf.defDim(ncid,'z',length(depth));
    dimid_time = netcdf.defDim(ncid,'time',length(time));
    
    % Create variable attributes
    varid_x = netcdf.defVar(ncid,'x','double',dimid_x);
    netcdf.putAtt(ncid,varid_x,'long_name','x');
    netcdf.putAtt(ncid,varid_x,'units','east_units');
    
    varid_y = netcdf.defVar(ncid,'y','double',dimid_y);
    netcdf.putAtt(ncid,varid_y,'long_name','y');
    netcdf.putAtt(ncid,varid_y,'units','north_units');
    
    varid_z = netcdf.defVar(ncid,'z','double',dimid_z);
    netcdf.putAtt(ncid,varid_z,'long_name','z');
    netcdf.putAtt(ncid,varid_z,'units','m');
    
    varid_time = netcdf.defVar(ncid,'time','double',dimid_time);
    netcdf.putAtt(ncid,varid_time,'long_name','time');
    netcdf.putAtt(ncid,varid_time,'units','days since 0');
    netcdf.putAtt(ncid,varid_time,'_CoordinateAxisType','Time');
    
    
    % variables
    varid_gamma_n = netcdf.defVar(ncid,'gamma_n','double',[dimid_x,dimid_y,dimid_z,dimid_time]);
    netcdf.putAtt(ncid,varid_gamma_n,'long_name','gamma_n');
    netcdf.putAtt(ncid,varid_gamma_n,'units','kgm3');
    
    netcdf.endDef(ncid)
    
    % Add data to file
    disp('Adding data to NetCDF file');
    netcdf.putVar(ncid,varid_x,x);
    netcdf.putVar(ncid,varid_y,y);
    netcdf.putVar(ncid,varid_z,depth);
    netcdf.putVar(ncid,varid_time,time);
    netcdf.putVar(ncid,varid_gamma_n,gamma_n);
    
    netcdf.close(ncid)
    
    
    disp(['Completed ' num2str(year(yy))]);
    toc
end % End of 'for each year (V20 file)' loop