function [power_results, prefs] = OrganizeDataAndRun_ANOVA(prefs)
%going to start with fast version, then need to add checks for exclusion & trial
%equality

%load csv
csvdata = textscan(fopen(prefs.csv_file), '%s%s%s%s', 'delimiter', ',');

%organize data
header = {csvdata{1}{1}, csvdata{2}{1}, csvdata{3}{1}, csvdata{4}{1}};
sub_list = csvdata{1}(2:end);
ans_list = csvdata{2}(2:end);
f1_levels_list = csvdata{3}(2:end);
f2_levels_list = csvdata{4}(2:end);
ans_as_nums = zeros(size(ans_list));
for c = 1:length(ans_list)
    ans_as_nums(c) = str2num(cell2mat(ans_list(c)));
end

%convert subject names to integers
sub_names = unique(csvdata{1}(2:end), 'stable');
sub_nums = zeros(length(csvdata{1}(2:end)), 1);
nPilotSubs = length(sub_names);
for sn = 1:nPilotSubs
    sub_nums(strcmp(sub_list, sub_names(sn)), 1) = sn;
end

%convert factor levels to integers
%f1
f1_names = unique(csvdata{3}(2:end));
f1_levels = zeros(length(csvdata{3}(2:end)), 1);
f1_num_levels = length(f1_names);
for cn = 1:f1_num_levels
    f1_levels(strcmp(f1_levels_list, f1_names(cn)), 1) = cn;
end

%f2
f2_names = unique(csvdata{4}(2:end));
f2_levels = zeros(length(csvdata{4}(2:end)), 1);
f2_num_levels = length(f2_names);
for cn = 1:f2_num_levels
    f2_levels(strcmp(f2_levels_list, f2_names(cn)), 1) = cn;
end


%subs, scores, and conditions as numbers
prefs.data = [sub_nums, ans_as_nums, f1_levels, f2_levels];

%column headers from CSV data
prefs.header = header;

%what kind of design (between, within, mixed)
if length(unique(f1_levels(sub_nums==1))) == 1 && length(unique(f2_levels(sub_nums==1))) == 1
    prefs.between_within_mixed = 1;
    %see if asking for equal condition allocation
    if strcmp(prefs.condition_allocation, 'even')
        prefs.condition_allocation = repmat(1/(f1_num_levels*f2_num_levels), 1, f1_num_levels*f2_num_levels);
    end
    

elseif length(unique(f1_levels(sub_nums==1))) == 1 || length(unique(f2_levels(sub_nums==1))) == 1
    
    prefs.between_within_mixed = 3;
    %make factor 1 the between subjects factor (if it isn't already)
    if length(unique(f1_levels(sub_nums==1))) > 1
        prefs.data = [sub_nums, ans_as_nums, f2_levels, f1_levels];
        tmp_names = f1_names;
        tmp_num_levels = f1_num_levels;
        f1_names = f2_names;
        f1_num_levels = f2_num_levels;
        f2_names = tmp_names;
        f2_num_levels = tmp_num_levels;
        prefs.header(3:4) = prefs.header([4,3]);
    end
    %see if asking for equal condition allocation
    if strcmp(prefs.condition_allocation, 'even')
        prefs.condition_allocation = repmat(1/(f1_num_levels), 1, f1_num_levels);
    end
    
elseif length(unique(f2_levels(sub_nums==1))) > 1 && length(unique(f2_levels(sub_nums==1))) > 1
    prefs.between_within_mixed = 2;
else
    error('Error in Data File')
end

%condition names (original as text)
prefs.f1_names = f1_names;
prefs.f2_names = f2_names;
prefs.f1_num_levels = f1_num_levels;
prefs.f2_num_levels = f2_num_levels;
prefs.sub_nums = sub_nums;
nConds = prefs.f1_num_levels * prefs.f2_num_levels;

%check to see whether each subject has an equal number of trials per
%condition (equal for condition, and for each subject)
%also, get subject means for each condition
sub_trial_counts = nan(nPilotSubs, nConds);
sub_means = nan(nPilotSubs, nConds);
for s = 1:nPilotSubs
    c = 0;
    for f1 = 1:prefs.f1_num_levels
        for f2 = 1:prefs.f2_num_levels
            c = c + 1;
            sub_trial_counts(s, c) = sum(prefs.data(:,1) == s & prefs.data(:,3) == f1 & prefs.data(:,4) == f2);
            sub_means(s,c) = mean(prefs.data(prefs.data(:,1) == s & prefs.data(:,3) == f1 & prefs.data(:,4) == f2,2));
        end
    end
end
nonZero = sub_trial_counts(sub_trial_counts > 0);
equal_trial_counts = all(nonZero == nonZero(1));
prefs.sub_trial_counts = sub_trial_counts;
prefs.sub_means = sub_means;

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
    %see if we can do fast sim (depends on subjects/trials/conditions)
    tMax = max(prefs.trial_range);
    nMax = max(prefs.N_range);

    if prefs.between_within_mixed == 2  %within-subjects design
        try
            biggest_array = zeros(nMax, (tMax+1)^nConds);
            if numel(biggest_array) <= prefs.max_array_size
                fast_sim = true;
            else
                fast_sim = false;
            end
        catch
            fast_sim = false;
        end
    elseif prefs.between_within_mixed == 3  %mixed-factor design
        try
            biggest_array = zeros(nMax, (tMax+1)^prefs.f2_num_levels);
            if numel(biggest_array) <= prefs.max_array_size
                fast_sim = true;
            else
                fast_sim = false;
            end
        catch
            fast_sim = false;
        end
    elseif prefs.between_within_mixed == 1 %between-subjects design
         try
            biggest_array = zeros(nMax, (tMax+1)^1);
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

% disp('Loading Data...');
% 
% %make plot of pilot data and simulation parameters
prefs = MakePilotGraph_ANOVA(prefs);

if prefs.show_pilot_data_only
    power_results = [];
    return
end

%main power analysis
if prefs.between_within_mixed == 2
    if fast_sim && equal_trial_counts && ~any_exclusion_criteria && ~varied_sim_trials
        power_results = PowerAnalysis_ANOVA_Within_Fast(prefs);
    elseif equal_trial_counts
        power_results = PowerAnalysis_ANOVA_Within_EqualTrials(prefs);
    else
        power_results = PowerAnalysis_ANOVA_Within_UnequalTrials(prefs);
    end
elseif prefs.between_within_mixed == 1
    if fast_sim && equal_trial_counts && ~any_exclusion_criteria && ~varied_sim_trials
        power_results = PowerAnalysis_ANOVA_Between_Fast(prefs);
    elseif equal_trial_counts
        power_results = PowerAnalysis_ANOVA_Between_EqualTrials(prefs);
    else
        power_results = PowerAnalysis_ANOVA_Between_UnequalTrials(prefs);
    end
elseif prefs.between_within_mixed == 3
    if fast_sim && equal_trial_counts && ~any_exclusion_criteria && ~varied_sim_trials
        power_results = PowerAnalysis_ANOVA_Mixed_Fast(prefs);
    elseif equal_trial_counts
        power_results = PowerAnalysis_ANOVA_Mixed_EqualTrials(prefs);
    else
        power_results = PowerAnalysis_ANOVA_Mixed_UnequalTrials(prefs);
    end
end

end


