function [regional_settings] = regional_settings_variance_plot(region)



%% Region-specific plot settings
if region == 1
    regional_settings.region_name = 'Region 1';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-42 -8];
    regional_settings.ylim = [52 74];
    
elseif region == 2
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
    
elseif region == 3
    regional_settings.region_name = 'Region 3';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-44 -20];
    regional_settings.ylim = [29 44];
    
elseif region == 4
    regional_settings.region_name = 'Region 4';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-63 -55];
    regional_settings.ylim = [40 45];
    
elseif region == 5
    regional_settings.region_name = 'Region 5';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-79 -43];
    regional_settings.ylim = [22 39];
    
elseif region == 6
    regional_settings.region_name = 'Region 6';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-40 -19];
    regional_settings.ylim = [11 21];
    
    
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
    
elseif region == 9
    regional_settings.region_name = 'Region 9';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-1 21];
    regional_settings.ylim = [-38 -19];
    
elseif region == 10
    regional_settings.region_name = 'Region 10';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-49 -40];
    regional_settings.ylim = [-31 -23];
    
elseif region == 11
    regional_settings.region_name = 'Region 11';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-42 -28];
    regional_settings.ylim = [-24 -18];
    
elseif region == 12
    regional_settings.region_name = 'Region 12';
    regional_settings.sal_caxis = ([10e-7 1]);
    regional_settings.temp_caxis = [10e-7 1];
    %regional_settings.sal_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 100];
    %regional_settings.temp_contours = [10e-12 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 100];
    regional_settings.sal_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e-0];
    regional_settings.temp_labels = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1 10e0];
    regional_settings.Scolours = (pmkmp(6,'CubicL'));
    regional_settings.Tcolours = (pmkmp(6,'CubicL'));
    regional_settings.xlim = [-56 -52];
    regional_settings.ylim = [-40 -36];
end



