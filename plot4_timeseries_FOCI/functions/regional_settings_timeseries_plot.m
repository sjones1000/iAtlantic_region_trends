function [regional_settings] = regional_settings_TS_plot(region,bathy)



%% Subset settings (How are we subdividing the water column for annual mean)


%% Which water masses to include on TS plots




%% Region-specific plot settings
if region == 1
    regional_settings.region_name = 'Subpolar MA Ridge Iceland';
    regional_settings.xlim = [-42 -8];
    regional_settings.ylim = [52 74];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [27.5 27.6; 27.625 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'ENACW','LSW','ISOW','DSOW'};
    regional_settings.boundary_polygon = [-41 -41 -32 -26 -14 -27; 53 60 65 63 61 53]; % 2 rows; x then y
    
elseif region == 2
    regional_settings.region_name = 'Rockall Trough';
    regional_settings.xlim = [-22 -3];
    regional_settings.ylim = [46 61];
    regional_settings.bathy_mask = (bathy > 1000);
    regional_settings.T_thresholds = [];
    regional_settings.S_thresholds = [];
    regional_settings.sig0_thresholds = []; 
    regional_settings.Yn_thresholds = [27.15 27.4; NaN NaN; 27.65 27.8; 27.9 35];
    regional_settings.WM_names = {'ENACW','MW','LSW','NEABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 3
    regional_settings.region_name = 'Central MA Ridge';
    regional_settings.xlim = [-44 -20];
    regional_settings.ylim = [29 44];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [26.9 27.2; NaN NaN; 27.65 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'ENACW','MW','uNADW','lNADW','NEABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 4
    regional_settings.region_name = 'Deep Sea Canyons_NW Atlantic';
    regional_settings.xlim = [-63 -55];
    regional_settings.ylim = [40 45];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [26.4 27.0; 27.0 27.5; 27.7 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'SBJW','ENACW','uNADW','lNADW','NEABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 5
    regional_settings.region_name = 'Sargasso Sea';
    regional_settings.xlim = [-79 -43];
    regional_settings.ylim = [22 39];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [25.8 26.6; 26.8 27.1; 27.1 27.4; NaN NaN; 27.7 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'STUW','WNACW','ENACW','AAIW','uNADW','lNADW','NEABW'};
    
elseif region == 6
    regional_settings.region_name = 'Eastern Tropical NA';
    regional_settings.xlim = [-40 -19];
    regional_settings.ylim = [11 21];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [NaN NaN; 24.6 25.4; 25.6 26.1; 26.1 26.6; 26.7 26.95; 27.05 27.4; 27.4 27.6; 27.9 35];
    regional_settings.WM_names = {'TSW','STUW','WNACW','ENACW','AAIW','uNADW','lNADW','NEABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 7
    regional_settings.region_name = 'Equatorial ocean fracture zones';
    regional_settings.xlim = [-26 -11.5];
    regional_settings.ylim = [-3 3.3];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [0 23.4; 23.4 25.9; 26.2 26.4; 26.4 26.6; 26.8 27.2; 27.5 27.7; 27.7 27.9; 27.9 35];
    regional_settings.WM_names = {'TSW','STUW','WSACW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 8
    regional_settings.region_name = 'Angola to Congo Lobe';
    regional_settings.xlim = [5 12.5];
    regional_settings.ylim = [-10.5 -4.5];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [22.3 25.8; 26.1 26.4; 26.4 26.7; 26.9 27.35; 27.5 27.7; 27.7 27.9; 27.9 35];
    regional_settings.WM_names = {'STUW','WSACW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 9
    regional_settings.region_name = 'Benguela Current';
    regional_settings.xlim = [-1 21];
    regional_settings.ylim = [-38 -19];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [26.1 26.4; 26.4 26.6; 26.8 27.3; 27.75 27.85; 27.85 27.95; 27.95 35];
    regional_settings.WM_names = {'WSACW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 10
    regional_settings.region_name = 'Brazil slope';
    regional_settings.xlim = [-49 -40];
    regional_settings.ylim = [-31 -23];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [24.7 25.4; 25.8 26.1; 26.1 26.5; 26.7 27.2; 27.65 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'STUW','WSACW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 11
    regional_settings.region_name = 'Vitoria-Trindade Seamount chain';
    regional_settings.xlim = [-42 -28];
    regional_settings.ylim = [-24 -18];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [24.8 25.6; 26 26.3; 26.3 26.5; 26.8 27.3; 27.6 27.75; 27.75 27.9; 27.9 35];
    regional_settings.WM_names = {'STUW','WSACW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
    
elseif region == 12
    regional_settings.region_name = 'Malvinas upwelling current_Argentina';
    regional_settings.xlim = [-56 -52];
    regional_settings.ylim = [-40 -36];
    regional_settings.bathy_mask = (bathy > 0);
    regional_settings.sig0_thresholds = [];
    regional_settings.Yn_thresholds = [25.7 26.5; 26.5 26.8; 26.9 27.3; 27.7 27.8; 27.8 27.9; 27.9 35];
    regional_settings.WM_names = {'SASW','ESACW','AAIW','uNADW','lNADW','AABW'};
    regional_settings.boundary_polygon = [];
end