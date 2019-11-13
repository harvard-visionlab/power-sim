function power_results = PowerAnalysis_ANOVA_Between_UnequalTrials(prefs)
if sum(prefs.condition_allocation) ~= 1
    error('Condition allocation does not add up to 1 (100%). Use fractions if this is due to rounding error (e.g., use 1/3 instead of .33)')
elseif length(prefs.condition_allocation) ~= prefs.f1_num_levels * prefs.f2_num_levels
    error('Must have an allocation percentage for each unique combination of factor levels.')
end

%simulation info
nSims = prefs.nSims; %number of experiments to simulate

cond = 0;
for f1 = 1:prefs.f1_num_levels
    for f2 = 1:prefs.f2_num_levels
        cond = cond + 1;
        level_f1(cond) = f1;
        level_f2(cond) = f2;
    end
end

maxTrialsNum = max(max(prefs.sub_trial_counts));

for f1 = 1:prefs.f1_num_levels
    for f2 = 1:prefs.f2_num_levels
        sub_nums{f1,f2} = unique(prefs.sub_nums(prefs.data(:,3) == f1 & prefs.data(:,4) == f2))';
        nPilotSubs(f1, f2) = length(sub_nums{f1,f2});
        pilot_data{f1,f2} = nan(nPilotSubs(f1,f2), maxTrialsNum);
        for s = 1:nPilotSubs(f1,f2)
            nTrialsInPilot{f1, f2}(s,1) = sum(prefs.data(:,1) == sub_nums{f1,f2}(s));
            pilot_data{f1,f2}(s,1:nTrialsInPilot{f1, f2}(s,1)) = prefs.data(prefs.data(:,1) == sub_nums{f1,f2}(s) ,2)';
        end
    end
end

sub_vector = prefs.N_range; %number of subs per simulation
prefs.trial_range = fliplr(prefs.trial_range);
trial_vector = prefs.trial_range(1,:); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
requested_sample_size = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = cell(length(trial_vector), length(sub_vector));
nConds = prefs.f1_num_levels*prefs.f2_num_levels;


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
        
        %number of subjects to simulate per condition
        nSubs_Total = sub_vector(sub_count);
        nSubs = zeros(1,nConds);
        
        for c = 1:nConds
            nSubs(c) = round(nSubs_Total*prefs.condition_allocation(c));
            cond_scores{c} = [];
        end
        
        subs_by_cond{trial_count, sub_count} = nSubs;
        requested_sample_size(trial_count, sub_count) = nSubs_Total;
        sample_size(trial_count, sub_count) = sum(nSubs);
        num_trials{trial_count, sub_count} = t;
        excluded_subs_total = 0;
        
        
        
        for c = 1:nConds
            total_num_trials = nSubs(c)*nSims*t(c);
            total_num_subs(c) = nSubs(c)*nSims;
            sim_ratio = total_num_trials/prefs.max_array_size;
            subs_per_round = floor(total_num_subs(c)/sim_ratio);
            [f1, f2] = Cond2FactorLevels(c, prefs);
            
            
            while size(cond_scores{c}, 1) < total_num_subs(c)
                
                %select random subjects
                sim_subs = randsample(1:nPilotSubs(f1,f2), subs_per_round, 'true');
                sim_subs = repelem(sim_subs, t(c))';
                
                cond_trials_per_sub = nTrialsInPilot{f1,f2}(sim_subs, 1);
                trial_nums = ceil(rand(length(sim_subs), 1).* cond_trials_per_sub);
                
                cond_scores_tmp = pilot_data{f1,f2}(sub2ind(size(pilot_data{f1,f2}), sim_subs, trial_nums));
                cond_scores_tmp2 = reshape(cond_scores_tmp, t(c), subs_per_round)';
                cond_scores_tmp3 = mean(cond_scores_tmp2, 2);
                
                included_subs = cond_scores_tmp3 > prefs.exclusion_min & cond_scores_tmp3 < prefs.exclusion_max;
                excluded_subs_total = excluded_subs_total + sum(~included_subs);
                cond_scores{c} = [cond_scores{c}; cond_scores_tmp3(included_subs)];
            end
            
            cond_scores{c} = cond_scores{c}(1:total_num_subs(c));
            cond_scores{c} = reshape(cond_scores{c}, nSubs(c), nSims);
            
        end
        
        % do condition comparisons
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            c1 = cond_scores{prefs.comps(comp,1)};
            c2 = cond_scores{prefs.comps(comp,2)};
            
            ds_vect{trial_count, sub_count}{comp} = (mean(c1) - mean(c2)) ./...
                (((nSubs(prefs.comps(comp,1)) - 1)*(std(c1).^2) + (nSubs(prefs.comps(comp,2)) - 1)*(std(c2).^2))/(nSubs(prefs.comps(comp,1)) + nSubs(prefs.comps(comp,2)) - 2)).^.5;
            [~,p] = ttest2(c1,c2);
            
            
            power_marker(comp, :) = p < prefs.alpha & ds_vect{trial_count, sub_count}{comp} > 0;
        end
        %%%
        
         if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %anova part
            Y = [];
            F1 = [];
            F2 = [];
           
            
            for c = 1:nConds
                [f1,f2] = Cond2FactorLevels(c, prefs);
                Y = [Y;cond_scores{c}];
                F1 = [F1; repmat(f1, nSubs(c), 1)];
                F2 = [F2; repmat(f2, nSubs(c), 1)];     
            end
                
            %between subjects anova
            p = anovan_matrix(Y, [F1, F2], 'model','interaction');
            
            if prefs.sig_ME1
                power_marker(end+1,:) = p(1,:) < prefs.alpha;
            end
            if prefs.sig_ME2
                power_marker(end+1,:) = p(2,:) < prefs.alpha;
            end
            if prefs.sig_int
                power_marker(end+1,:) = p(3,:) < prefs.alpha;
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
power_results.requested_n = requested_sample_size; %requested sample size for each simulated design
power_results.num_trials = num_trials; %number of trials for each design
%power_results.dz_vect = dz_vect; %effect size vector for each design
power_results.sub_vector = sample_size(1,:);
power_results.trial_vector = trial_vector;
power_results.excluded_subs_avg = excluded_subs;
power_results.subs_by_cond = subs_by_cond;  %subs in each cond for each design

end