function [regional_settings] = regional_settings_variance_plot(region)



%% Region-specific plot settings

if region == 2
    regional_settings.region_name = 'Rockall Trough';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-22 -3];
    regional_settings.ylim = [46 61];
    
elseif region == 7
    regional_settings.region_name = 'Region 7';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-26 -11.5];
    regional_settings.ylim = [-3 3.3];
    
elseif region == 8
    regional_settings.region_name = 'Region 8';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [5 12.5];
    regional_settings.ylim = [-10.5 -4.5];
end



