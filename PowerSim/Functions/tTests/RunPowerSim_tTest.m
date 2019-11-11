function power_results = RunPowerSim_tTest(prefs)
%main function for t-test power analyses

%load additional settings
prefs = MainSetup(prefs);

%load data, organize data, run power analysis
[power_results, prefs] = OrganizeDataAndRun_tTests(prefs);

%exit if only showing pilot data
if prefs.show_pilot_data_only
    warning('Show pilot data only mode: Exiting before simulation.');
    return
end


%plot the heat map
PlotHeatMap(power_results, prefs);

end