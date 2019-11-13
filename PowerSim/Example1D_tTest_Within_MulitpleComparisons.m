clear
close all
addpath(genpath('Functions'))

%path to csv file containing trial-level data
prefs.csv_file = 'Example_Data/Data_tTest_Within_Multi.csv';

%show pilot data only (set to true if you want to view pilot data but not
%simulate design...useful for making sure you set everything up right
%must be set to false to run actual power analysis
prefs.show_pilot_data_only = false;

%interval of N to simulate (e.g, 10-100 by 10)
%for between-subjects designs, this is TOTAL subjects (not per condition)
prefs.N_range = 50:50:250;

%interval of trials per condition to simulate (e.g, 8-24 by 4)
%can either be a single row (where trial amount applied to every condition)
prefs.trial_range = 8:8:24;
%or variable trials per condition (e.g., 2 conditions with 2:1 trial ratio)
%need a row for each condition, in this case
%prefs.trial_range = [8:4:24; 4:2:12];

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
prefs.comps = [4, 1
    4,2
    4,3
    2,1
    3,1];

%FOR BETWEEN-SUBJECTS DESIGNS ONLY (code ignores value otherwise)
%how participants should be split between conditions, must sum to 1 (100%)
%for example, if 60 participants in 2 condition between-subjects design,
%prefs.condition_allocation = [.5, .5] would have 30 subs/condition.
%[.75, .25] would result in condition 1 = 45 subs, condition 2 = 15 subs
%[1/3, 2/3] would result in condition 1 = 20 subs, condition 2 = 40 subs
%if you want subjects to be split evenly between all your conditions, you
%can assign: prefs.condition_allocation = 'even';
prefs.condition_allocation = 'even';

%subject-level exlcusion info
%if DV is percent correct, use decimal (e.g., .6)
prefs.exclusion_min = []; 
prefs.exclusion_max = [];

%Run Power Analysis with these settings
power_results = RunPowerSim_tTest(prefs);
