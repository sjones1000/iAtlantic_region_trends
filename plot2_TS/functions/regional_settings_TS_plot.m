function [regional_settings] = regional_settings_TS_plot(region,bathy)



%% Subset settings (How are we subdividing the water column for annual mean)

% Mask data using bathymetry as a threshold?
if region == 2
regional_settings.bathy_mask = (bathy > 1000);
elseif region == 7
    regional_settings.bathy_mask = (bathy > 0);
elseif region == 8
    regional_settings.bathy_mask = (bathy > 0);
end


%% Region-specific plot settings

if region == 2
    regional_settings.region_name = 'Rockall Trough';
elseif region == 7
    regional_settings.region_name = 'Region 7';
elseif region == 8
    regional_settings.region_name = 'Region 8';
end
