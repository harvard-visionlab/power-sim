function power_results = PowerAnalysis_ANOVA_Within_EqualTrials(prefs)

%simulation info
nSims = prefs.nSims; %number of experiments to simulate
nPilotSubs = length(unique(prefs.data(:,1))); %how many subjects in actual data
sub_vector = prefs.N_range; %number of subs per simulation
prefs.trial_range = fliplr(prefs.trial_range);
trial_vector = prefs.trial_range(1,:); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest
nConds = prefs.f1_num_levels*prefs.f2_num_levels; %number of conditions
nPilotTrials = sum(prefs.data(:,1) == 1 & prefs.data(:,3) == 1 & prefs.data(:,4) == 1);

%organize data sub*trial*cond
pilot_data = nan(nPilotSubs, nPilotTrials, nConds);
for s = 1:nPilotSubs
    cond = 0;
    for f1 = 1:prefs.f1_num_levels
        for f2 = 1:prefs.f2_num_levels
            cond = cond + 1;
            pilot_data(s,:,cond) = prefs.data(prefs.data(:,1) == s & prefs.data(:,3) == f1 & prefs.data(:,4) == f2, 2)';
        end
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
    
    %calculate measurement variability for each subject based upon
    %condition means and number of trials
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        pc = round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)));
        disp([num2str(pc), '% Complete']);
        
        
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
        
        for c = 1:nConds
            cond_scores{c} = sub_scores(:,:,c);
        end
        
        % do condition comparisons
        %p = cell(1, nComps);
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            diff_scores = cond_scores{prefs.comps(comp,1)} - cond_scores{prefs.comps(comp,2)};
            dz_vect{trial_count, sub_count}{comp} = mean(diff_scores)./std(diff_scores);
            [~,p] = ttest(diff_scores);
            power_marker(comp, :) = p < prefs.alpha & dz_vect{trial_count, sub_count}{comp} > 0;
        end
        
        if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %anova part
            
            Y = [];
            F1 = [];
            F2 = [];
            c = 0;
            for f1 = 1:prefs.f1_num_levels
                F1 = [F1; f1*ones(nSubs * prefs.f2_num_levels,1)];
                for f2 = 1:prefs.f2_num_levels
                    c = c + 1;
                    Y = [Y; cond_scores{c}];
                    F2 = [F2; f2*ones(nSubs,1)];
                end
            end
            
            S = repmat(1:nSubs, 1, nConds)';
            
            stats = rm_anova2_matrix(Y,S,F1,F2);
            
            if prefs.sig_ME1
                power_marker(end+1,:) = stats.pA < prefs.alpha;
            end
            if prefs.sig_ME2
                power_marker(end+1,:) = stats.pB < prefs.alpha;
            end
            if prefs.sig_int
                power_marker(end+1,:) = stats.pAB < prefs.alpha;
            end
        end
        
        power(trial_count, sub_count) = mean(all(power_marker, 1));
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