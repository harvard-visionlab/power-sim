function [power_results, prefs] = OrganizeDataAndRun_tTests(prefs)

%load csv
csvdata = textscan(fopen(prefs.csv_file), '%s%s%s', 'delimiter', ',');

%organize data
header = {csvdata{1}{1}, csvdata{2}{1}, csvdata{3}{1}};
sub_list = csvdata{1}(2:end);
ans_list = csvdata{2}(2:end);
cond_list = csvdata{3}(2:end);
ans_as_nums = str2num(cell2mat(ans_list));

%convert subject names to integers
sub_names = unique(csvdata{1}(2:end), 'stable');
sub_nums = zeros(length(csvdata{1}(2:end)), 1);
for sn = 1:length(sub_names)
    sub_nums(strcmp(sub_list, sub_names(sn)), 1) = sn;
end

%convert condition names to integers
cond_names = unique(csvdata{3}(2:end));
cond_nums = zeros(length(csvdata{3}(2:end)), 1);
for cn = 1:length(cond_names)
    cond_nums(strcmp(cond_list, cond_names(cn)), 1) = cn;
end

%see if asking for equal condition allocation
if strcmp(prefs.condition_allocation, 'even')
   prefs.condition_allocation = repmat(1/length(cond_names), 1, length(cond_names));
end

%subs, scores, and conditions as numbers
prefs.data = [sub_nums, ans_as_nums, cond_nums];

%condition names (original as text)
prefs.cond_names = cond_names;

%column headers from CSV data
prefs.header = header;

%determine whether a within-subject or between-subject design
subs = unique(prefs.data(:,1));
if length(unique(prefs.data(prefs.data(:,1) == subs(1),3))) == 1
    prefs.within_between = 2; %between subjects design
else
    prefs.within_between = 1; %within subjects design
end

%check to see whether each subject has an equal number of trials per
%condition (equal for condition, and for each subject)
nConds = length(cond_names);
sub_trial_counts = nan(length(subs), nConds);
for s = 1:length(subs)
    for c = 1:nConds
        sub_trial_counts(s, c) = sum(prefs.data(:,1) == s & prefs.data(:,3) == c);
    end
end
nonZero = sub_trial_counts(sub_trial_counts > 0);
equal_trial_counts = all(nonZero == nonZero(1));
prefs.sub_trial_counts = sub_trial_counts;

%check to see if variable amount of trials per condition
if size(prefs.trial_range, 1) == 1
    varied_sim_trials = false;
elseif size(prefs.trial_range, 1) > 1 && size(prefs.trial_range, 1) == nConds
    varied_sim_trials = true;
else
    error('Error with trial range. Must be single row of values (to be applied to all conditions), or a matrix with with a row for each condition');
end
prefs.varied_sim_trials = varied_sim_trials;

%determine whether DV is percent correct (only 2 choices)
if length(unique(prefs.data(:,2))) == 2 && ~varied_sim_trials
    
    tMax = max(prefs.trial_range);
    nMax = max(prefs.N_range);
    
    if prefs.within_between == 1  %within-subjects design
        
        try
           biggest_array = zeros(nMax, (tMax+1)^nConds); %see if matlab can handle data of this size
           if numel(biggest_array) <= prefs.max_array_size
                fast_sim = true;
            else
                fast_sim = false;
            end
        catch
            fast_sim = false;
        end
    elseif prefs.within_between == 2  %between-subjects design design
        try
            biggest_array = zeros(nMax, (tMax+1)^1); %see if matlab can handle data of this size (should be fine for between-subs design)
            if numel(biggest_array) <= prefs.max_array_size
                fast_sim = true;
            else
                fast_sim = false;
            end
        catch
            fast_sim = false;
        end
    end
else
    fast_sim = false;
end

any_exclusion_criteria = ~isempty(prefs.exclusion_min) || ~isempty(prefs.exclusion_max);
if isempty(prefs.exclusion_min)
    prefs.exclusion_min = -inf;
end
if isempty(prefs.exclusion_max)
    prefs.exclusion_max = inf;
end

%make pilot data graph
MakePilotGraph_tTests(prefs);

if prefs.show_pilot_data_only
    power_results = [];
    return
end

%main power analysis
if prefs.within_between == 1
    if fast_sim && equal_trial_counts && ~any_exclusion_criteria && ~varied_sim_trials
        power_results = PowerAnalysis_tTest_Within_Fast(prefs);
    elseif equal_trial_counts
        power_results = PowerAnalysis_tTest_Within_EqualTrials(prefs);
    else
        power_results = PowerAnalysis_tTest_Within_UnequalTrials(prefs);
    end
elseif prefs.within_between == 2
    if fast_sim && equal_trial_counts && ~any_exclusion_criteria && ~varied_sim_trials
        power_results = PowerAnalysis_tTest_Between_Fast(prefs);
    elseif equal_trial_counts
        power_results = PowerAnalysis_tTest_Between_EqualTrials(prefs);
    else
        power_results = PowerAnalysis_tTest_Between_UnequalTrials(prefs);
    end
end

end


