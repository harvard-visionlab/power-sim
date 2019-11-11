function power_results = PowerAnalysis_tTest_Within_EqualTrials(prefs)

%simulation info
nSims = prefs.nSims; %number of experiments to simulate
nPilotSubs = length(unique(prefs.data(:,1))); %how many subjects in actual data
sub_names = unique(prefs.data(:,1));
sub_vector = prefs.N_range; %number of subs per simulation
prefs.trial_range = fliplr(prefs.trial_range);
trial_vector = prefs.trial_range(1,:); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest
nConds = length(unique(prefs.data(:,3))); %number of conditions
cond_names = unique(prefs.data(:,3));
nPilotTrials = sum(prefs.data(:,1) == sub_names(1) & prefs.data(:,3) == cond_names(1));

%organize data sub*trial*cond
pilot_data = zeros(nPilotSubs, nPilotTrials, nConds);
for s = 1:nPilotSubs
    for cond = 1:nConds
        pilot_data(s,:,cond) = prefs.data(prefs.data(:,1) == s & prefs.data(:,3) == cond, 2)';
    end
end

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = cell(length(trial_vector), length(sub_vector));
dz_vect = cell(1, nComps);

for trial_count= 1:length(trial_vector)
    
    if prefs.varied_sim_trials
        t = prefs.trial_range(:, trial_count);
    else
        t = repmat(trial_vector(trial_count), nConds, 1);
    end
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        disp([num2str(round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)))), '% Complete']);
        
        %number of subjects to simulate
        nSubs = sub_vector(sub_count);
        sample_size(trial_count, sub_count) = nSubs;
        num_trials{trial_count, sub_count} = t;
        
        total_num_trials = nSubs*nSims*max(t);
        total_num_subs = nSubs*nSims;
        sim_ratio = total_num_trials/prefs.max_array_size;
        subs_per_round = floor(total_num_subs/sim_ratio);
        sub_scores = [];
        excluded_subs_total = 0;
        
        
        while size(sub_scores, 1) < total_num_subs
            
            %select random subjects
            sim_subs_total = randsample(1:nPilotSubs, subs_per_round, 'true');
            
            
            %select random trials, for each condition, generate data
            cond_scores_tmp3 = [];
            for c = 1:nConds
                sim_subs = repelem(sim_subs_total, t(c))';
                trial_nums = ceil(rand(length(sim_subs), 1)*nPilotTrials);
                cond_scores_tmp = pilot_data(sub2ind(size(pilot_data), sim_subs, trial_nums, repmat(c, length(sim_subs), 1)));

                cond_scores_tmp2 = reshape(cond_scores_tmp, t(c), subs_per_round)';
                cond_scores_tmp3(:,c) = mean(cond_scores_tmp2, 2);
            end
            
            
            t_ratio = t'/sum(t);
            weighted_average = sum(t_ratio .* cond_scores_tmp3,2);
            included_subs = weighted_average > prefs.exclusion_min & weighted_average < prefs.exclusion_max;
            excluded_subs_total = excluded_subs_total + sum(~included_subs);
            
            
            sub_scores = [sub_scores; cond_scores_tmp3(included_subs, :)];
            
            
        end
        
        sub_scores = sub_scores(1:total_num_subs, :);
        sub_scores = reshape(sub_scores, nSubs, nSims, nConds);
        
        % do condition comparisons
        %p = cell(1, nComps);
        power_marker{trial_count, sub_count} = zeros(nComps, nSims);
        for comp = 1:nComps
            diff_scores = sub_scores(:,:,prefs.comps(comp,1)) - sub_scores(:,:,prefs.comps(comp,2));
            dz_vect{trial_count, sub_count}{comp} = mean(diff_scores)./std(diff_scores);
            [~,p] = ttest(diff_scores);
            power_marker{trial_count, sub_count}(comp, :) = p < prefs.alpha & dz_vect{trial_count, sub_count}{comp} > 0;
        end
        
        power(trial_count, sub_count) = mean(all(power_marker{trial_count, sub_count}, 1));
        excluded_subs(trial_count, sub_count) = excluded_subs_total/nSims;
    end
end

clc
disp('100% Complete')

%output information
power_results.power = power; %power for each simulated design
power_results.n = sample_size; %sample size for each simulated design
power_results.num_trials = num_trials; %number of trials for each design
power_results.dz_vect = dz_vect; %effect size vector for each design
power_results.sub_vector = sub_vector;
power_results.trial_vector = trial_vector;
power_results.excluded_subs_avg = excluded_subs;


end