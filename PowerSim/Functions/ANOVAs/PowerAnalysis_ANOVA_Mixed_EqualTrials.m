function power_results = PowerAnalysis_ANOVA_Mixed_EqualTrials(prefs)
if sum(prefs.condition_allocation) ~= 1
    error('Condition allocation does not add up to 1 (100%). Use fractions if this is due to rounding error (e.g., use 1/3 instead of .33)')
elseif length(prefs.condition_allocation) ~= prefs.f1_num_levels
    error('Must have a condition allocation (prefs.condition_allocation) amount for each between-subjects level.')
end

%simulation info
nSims = prefs.nSims; %number of experiments to simulate
nPilotTrials = sum(prefs.data(:,1) == 1 & prefs.data(:,3) == 1 & prefs.data(:,4) == 1);

%number of subs for each between-subs factor level

num_between_levels = prefs.f1_num_levels;
num_within_levels =  prefs.f2_num_levels;

cond = 0;
for f1 = 1:num_between_levels
    for f2 = 1:num_within_levels
        cond = cond + 1;
        level_between(cond) = f1;
        level_within(cond) = f2;
    end
end


for c = 1:num_between_levels
    between_subs{c} = unique(prefs.sub_nums(prefs.data(:,3) == c))';
    nPilotSubs(c) = length(between_subs{c});
end

sub_vector = prefs.N_range; %number of subs per simulation
prefs.trial_range = fliplr(prefs.trial_range);
trial_vector = prefs.trial_range(1,:); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
requested_sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = cell(length(trial_vector), length(sub_vector));

%set up data to be easily sampled in simulations
%%%%
%organize data {between_cond}sub*trial*within_cond
pilot_data = cell(1, num_between_levels);
for f1 = 1:prefs.f1_num_levels
    pilot_data{f1} = nan(nPilotSubs(f1), nPilotTrials, prefs.f2_num_levels);
    for s = 1:nPilotSubs(f1)
        for f2 = 1:prefs.f2_num_levels
            pilot_data{f1}(s,:,f2) = prefs.data(prefs.data(:,1) == between_subs{f1}(s) & prefs.data(:,3) == f1 & prefs.data(:,4) == f2, 2)';
        end
    end
end


for trial_count= 1:length(trial_vector)
    
    if prefs.varied_sim_trials
        t = prefs.trial_range(:, trial_count);
    else
        t = repmat(trial_vector(trial_count), prefs.f1_num_levels*prefs.f2_num_levels, 1);
    end
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        pc = round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)));
        disp([num2str(pc), '% Complete']);
        
        nSubs_Total = sub_vector(sub_count);
        excluded_subs_total = 0;
        num_trials{trial_count, sub_count} = t;
        
        for f1 = 1:prefs.f1_num_levels
            %number of subjects to simulate
            
            nSubs(f1) = round(nSubs_Total*prefs.condition_allocation(f1));
            total_num_trials = nSubs(f1)*nSims*max(t);
            total_num_subs(f1) = nSubs(f1)*nSims;
            sim_ratio = total_num_trials/prefs.max_array_size;
            subs_per_round = floor(total_num_subs(f1)/sim_ratio);
            sub_scores = [];
            
            while size(sub_scores, 1) < total_num_subs(f1)
                
                %select random subjects
                sim_subs_total = randsample(1:nPilotSubs(f1), subs_per_round, 'true');
                
                
                %select random trials, for each condition, generate data
                cond_scores_tmp3 = [];
                for f2 = 1:prefs.f2_num_levels
                    c = FactorLevels2Cond(f1, f2, prefs);
                    sim_subs = repelem(sim_subs_total, t(c))';
                    trial_nums = ceil(rand(length(sim_subs), 1)*nPilotTrials);
                    
                    cond_scores_tmp = pilot_data{f1}(sub2ind(size(pilot_data{f1}), sim_subs, trial_nums, repmat(f2, length(sim_subs), 1)));
                    cond_scores_tmp2 = reshape(cond_scores_tmp, t(c), subs_per_round)';
                    cond_scores_tmp3(:,f2) = mean(cond_scores_tmp2, 2);
                end
                
                %which t values to use for this factor
                ti = FactorLevels2Cond(f1, 1, prefs);
                tf = FactorLevels2Cond(f1, prefs.f2_num_levels, prefs);
                t_now = t(ti:tf);
                
                t_ratio = t_now'/sum(t_now);
                weighted_average = sum(t_ratio .* cond_scores_tmp3,2);
                included_subs = weighted_average > prefs.exclusion_min & weighted_average < prefs.exclusion_max;
                excluded_subs_total = excluded_subs_total + sum(~included_subs);
                
                sub_scores = [sub_scores; cond_scores_tmp3(included_subs, :)];
            end
            
            sub_scores = sub_scores(1:total_num_subs(f1), :);
            sub_scores = reshape(sub_scores, nSubs(f1), nSims, prefs.f2_num_levels);
            
            for f2 = 1:prefs.f2_num_levels
                cond_scores{f1}{f2} = sub_scores(:,:,f2);
            end
        end
        
        sample_size(trial_count, sub_count) = sum(nSubs);
        subs_by_cond{trial_count, sub_count} = repelem(nSubs, prefs.f2_num_levels);
        requested_sample_size(trial_count, sub_count) = nSubs_Total;
        
        % do condition comparisons
        %p = cell(1, nComps);
        %%%
        % do t-tests of interest
        %%%
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            %figure out if a within or between comparison
            if level_between(prefs.comps(comp,1)) == level_between(prefs.comps(comp,2)) %within comparison
                bl = level_between(prefs.comps(comp,1));
                wl1 = level_within(prefs.comps(comp,1));
                wl2 = level_within(prefs.comps(comp,2));
                diff_scores = cond_scores{bl}{wl1} - cond_scores{bl}{wl2};
                diff_scores = reshape(diff_scores, nSubs(bl), nSims);
                dz_vect = mean(diff_scores)./std(diff_scores);
                [~,p] = ttest(diff_scores);
                power_marker(comp, :) = p < prefs.alpha & dz_vect > 0;
                
            else %between comparison
                bl1 = level_between(prefs.comps(comp,1));
                bl2 = level_between(prefs.comps(comp,2));
                wl1 = level_within(prefs.comps(comp,1));
                wl2 = level_within(prefs.comps(comp,2));
                
                c1 = cond_scores{bl1}{wl1};
                c2 = cond_scores{bl2}{wl2};
                ds_vect = (mean(c1) - mean(c2)) ./...
                    (((nSubs(bl1) - 1)*(std(c1).^2) + (nSubs(bl2) - 1)*(std(c2).^2))/(nSubs(bl1) + nSubs(bl2) - 2)).^.5;
                [~,p] = ttest2(c1,c2);
                power_marker(comp, :) = p < prefs.alpha & ds_vect > 0;
                
            end
        end
        %%%
        
        if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %anova part
            Y = [];
            S = [];
            WF = [];
            BF = [];
            sm = 1;
            
            for between_level = 1:num_between_levels
                
                for c = 1:num_within_levels
                    Y = [Y;reshape(cond_scores{between_level}{c}, nSubs(between_level), nSims)];
                    S = [S; (sm:sm+nSubs(between_level)-1)'];
                    WF = [WF; repmat(c, nSubs(between_level), 1)];
                    BF = [BF; repmat(between_level, nSubs(between_level), 1)];
                end
                sm = sm + nSubs(between_level);
            end
            
            stats = mixed_anova_matrix(Y,S,WF,BF);
            
            if prefs.sig_ME1
                power_marker(end+1,:) = stats.Pbs < prefs.alpha;
            end
            if prefs.sig_ME2
                power_marker(end+1,:) = stats.Pws < prefs.alpha;
            end
            if prefs.sig_int
                power_marker(end+1,:) = stats.Pint < prefs.alpha;
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
powre_results.subs_by_cond = subs_by_cond;

end