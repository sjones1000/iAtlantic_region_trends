function [regional_settings] = regional_settings_TS_plot(region,bathy)



%% Subset settings (How are we subdividing the water column for annual mean)


%% Which water masses to include on TS plots
%["ENACW" = 1;"WNACW" = 2;"ESACW" = 3;"WSACW" =4;"AAIW" =5;"SAIW" = 6;"MW" = 7;
%    "NADW" = 8;"LSW" = 9;"ISOW" = 10;"DSOW" = 11;"AABW" = 12;"CDW" = 13;"WSBW" = 14;"NEABW" = 15]




%% Region-specific plot settings
if region == 1
    regional_settings.region_name = 'Subpolar MA Ridge Iceland';
    regional_settings.xlim = [-42 -8];
    regional_settings.ylim = [52 74];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [-41 -41 -32 -26 -14 -27; 53 60 65 63 61 53]; % 2 rows; x then y
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 2
    regional_settings.region_name = 'Rockall Trough';
    regional_settings.xlim = [-22 -3];
    regional_settings.ylim = [46 61];
    regional_settings.bathy_mask = (bathy > 1000);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 3
    regional_settings.region_name = 'Central MA Ridge';
    regional_settings.xlim = [-44 -20];
    regional_settings.ylim = [29 44];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 4
    regional_settings.region_name = 'Deep Sea Canyons_NW Atlantic';
    regional_settings.xlim = [-63 -55];
    regional_settings.ylim = [40 45];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = [34.6 36 1.2 15.2]; % [x1 x2 y1 y2];
    
elseif region == 5
    regional_settings.region_name = 'Sargasso Sea';
    regional_settings.xlim = [-79 -43];
    regional_settings.ylim = [22 39];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 6
    regional_settings.region_name = 'Eastern Tropical NA';
    regional_settings.xlim = [-40 -19];
    regional_settings.ylim = [11 21];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 7
    regional_settings.region_name = 'Equatorial ocean fracture zones';
    regional_settings.xlim = [-26 -11.5];
    regional_settings.ylim = [-3 3.3];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 8
    regional_settings.region_name = 'Angola to Congo Lobe';
    regional_settings.xlim = [5 12.5];
    regional_settings.ylim = [-10.5 -4.5];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 9
    regional_settings.region_name = 'Benguela Current';
    regional_settings.xlim = [-1 21];
    regional_settings.ylim = [-38 -19];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 10
    regional_settings.region_name = 'Brazil slope';
    regional_settings.xlim = [-49 -40];
    regional_settings.ylim = [-31 -23];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 11
    regional_settings.region_name = 'Vitoria-Trindade Seamount chain';
    regional_settings.xlim = [-42 -28];
    regional_settings.ylim = [-24 -18];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
    
elseif region == 12
    regional_settings.region_name = 'Malvinas upwelling current_Argentina';
    regional_settings.xlim = [-56 -52];
    regional_settings.ylim = [-40 -36];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.water_mass_plot = [1:15];
    regional_settings.boundary_polygon = [];
    regional_settings.axis_limits = []; % [x1 x2 y1 y2];
end