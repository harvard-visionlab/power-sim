clear
close all
addpath(genpath('Functions'))

%path to csv file containing trial-level data
prefs.csv_file = 'Example_Data/Data_ANOVA_Between.csv';

%show pilot data only (set to true if you want to view pilot data but not
%simulate design...useful for making sure you set everything up right
%must be set to false to run actual power analysis
prefs.show_pilot_data_only = false;

%interval of N to simulate (e.g, 50-300 by 50)
prefs.N_range = 50:50:300; 

%interval of trials per condition to simulate (e.g, 8-20 by 4)
prefs.trial_range = 16:16:48;
%%say you want to simulate having a 1:2:1:2 ratio of 4 conditions
%prefs.trial_range = [12:12:24; 24:24:48; 12:12:24; 24:24:48];

%p value to use in statisical test during simulation
prefs.alpha = .05;

%number of experiments to simulate per trial*N combination
%higher number of sims will give more stable/accurate power estimates, 
%but will be slower. 10000 or 100000 is usually good
prefs.nSims = 10000;

%what comparisons do you want to make? Should be a comparison * 2 vector,
%with condition that should be larger on the left
%for example, if you expect condition 2 to be larger than condition 1, you
%should enter [2, 1];
%when you run this script, a graph will display how your conditions have
%been numbered (adjust below and run again if necessary)
prefs.comps = [2, 1];

%does the first main effect need to be significant to be a "success"
%note that for mixed-factors designs, between-subjects factor will always 
%be considered "factor 1"
prefs.sig_ME1 = true;

%does the second main effect need to be significant to be a "success"
%note that for mixed-factors designs, within-subjects factor will always 
%be considered "factor 2"
prefs.sig_ME2 = false;

%does the interaction need to be significant to be a "success"
prefs.sig_int = true;

%FOR BETWEEN-SUBJECTS OR MIXED DESIGNS ONLY (ignored otherwise)
%how participants should be split between between-factor levels
%needs a value for each between-subjects factor level, and sum to 1 (100%)
%For mixed designs, should have a value for each between-subjects factor
%level
%For between subjects designs, need a value for each condition (total
%number of conditions is #factor1levels*#factor2levels), e.g., for a 3x2
%between subjects design, would need 6 values.
%for example, if 60 participants in 2 condition between-subjects design,
%prefs.condition_allocation = [.5, .5] would have 30 subs/condition.
%[.75, .25] would result in condition 1 = 45 subs, condition 2 = 15 subs
%[1/3, 2/3] would result in condition 1 = 20 subs, condition 2 = 40 subs
prefs.condition_allocation = 'even';

%subject-level exlcusion criteria
prefs.exclusion_min = [.6]; %if DV is percent correct, use decimal (e.g., .6)
prefs.exclusion_max = [];

%Run Power Analysis with these settings
power_results = RunPowerSim_ANOVA(prefs);