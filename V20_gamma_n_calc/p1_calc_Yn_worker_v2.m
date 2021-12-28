function p1_calc_Yn_worker_v2(path,file,out_dir)
addpath(genpath('/gxfs_work1/geomar/smomw294/iAtlantic/Subprojects/iAtlantic_region_trends/other_functions'));

disp(['************** Working on Region ' file ' *******************']);
disp(['******************************************************']);
tic

% %% Read global nc vars
filename = [path,file];
z2d = ncread(filename,'gdept_0'); % partial grid cells at bottom due to bathymetry
x = double(ncread(filename,'x'));
y = double(ncread(filename,'y'));
lon2d = ncread(filename,'glamt');
lat2d = ncread(filename,'gphit');
dx2d = ncread(filename,'e1t');
dy2d = ncread(filename,'e2t');
bathyind = ncread(filename,'mbathy');
time = ncread(filename,'time');
depth = ncread(filename,'depth');

temp = ncread(filename,'votemper');
sal = ncread(filename,'vosaline');
tmask = ncread(filename,'tmask');

[xN,yN,zN,tN]=size(z2d);

tmask = repmat(tmask,1,1,1);
sal(tmask==0)=NaN; temp(tmask==0)=NaN;
clear tmask

i_min = 351; i_max = xN;
j_min = 201; j_max = yN; 

sal = sal(i_min:i_max,j_min:j_max,:,:);
temp = temp(i_min:i_max,j_min:j_max,:,:);

z2d = z2d(i_min:i_max,j_min:j_max,:).*-1;


% Need a 4d depth var
% 
lon4d = repmat(lon2d(i_min:i_max,j_min:j_max),1,1,zN);
lat4d = repmat(lat2d(i_min:i_max,j_min:j_max),1,1,zN);
x=x(i_min:i_max);
y=y(j_min:j_max);
% Calc pressure
% p = gsw_p_from_z(z,lat);
press = gsw_p_from_z(z2d,lat4d);
gamma_n = temp.*NaN;

% parpool(4)
% maxNumCompThreads(4)
% parfor (mm=1:tN,4)
for mm=1:tN
	% Calc annual mean
        T3d = squeeze(temp(:,:,:,mm));
        S3d = squeeze(sal(:,:,:,mm));
        z3d = squeeze(press(:,:,:));
        y3d = squeeze(lat4d(:,:,:));
        x3d = squeeze(lon4d(:,:,:));

	gamma_n(:,:,:,mm) = eos80_legacy_gamma_n(S3d,T3d,...
        z3d,x3d,y3d);
        disp(['Month ' num2str(mm) ])
end
    %% Write netcdf
    fileout = [out_dir,file(1:end-3),'_4.nc']; % <- Save data location
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
    
%     poolobj = gcp('nocreate')
%     delete(poolobj);
    
    disp(['Completed ' file]);
    toc
end % End of 'for each year (V20 file)' loop
