function power_results = RunPowerSim_ANOVA(prefs)
%main function for ANOVA power analyses

%load additional settings
prefs = MainSetup(prefs);

%load data, organize data, run power analysis
[power_results, prefs] = OrganizeDataAndRun_ANOVA(prefs);

if prefs.show_pilot_data_only
    warning('Show pilot data only mode: Exiting before simulation.');
    return
end

PlotHeatMap(power_results, prefs);

end