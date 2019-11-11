function power_results = PowerAnalysis_ANOVA_Mixed_Fast(prefs)

if sum(prefs.condition_allocation) ~= 1
    error('Condition allocation does not add up to 1 (100%). Use fractions if this is due to rounding error (e.g., use 1/3 instead of .33)')
elseif length(prefs.condition_allocation) ~= prefs.f1_num_levels
    error('Must have an allocation amount for each between-subjects level.')
end

%simulation info
nSims = prefs.nSims; %number of experiments to simulate

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
trial_vector = fliplr(prefs.trial_range); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
requested_sample_size = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = zeros(length(trial_vector), length(sub_vector));

%get subject-level long-form data for simulation
prefs.subs = [];
prefs.F1 = [];
prefs.F2 = [];
prefs.Y = [];

for f1 = 1:num_between_levels
    for sub = 1:nPilotSubs(f1)
        for f2 = 1:num_between_levels
            prefs.subs = [prefs.subs; between_subs{f1}(sub)];
            prefs.F1 = [prefs.F1; f1];
            prefs.F2 = [prefs.F2; f2];
            prefs.Y = [prefs.Y; mean(prefs.data(prefs.data(:,1) == between_subs{f1}(sub) & prefs.data(:,4) == f2,2))];
        end
    end
end


for trial_count= 1:length(trial_vector)
    
    t = trial_vector(trial_count);
    
    
    for between_level = 1:num_between_levels
        %determine condition difference pdf for each subject in pilot data
        cFinal{between_level} = zeros(nPilotSubs(between_level), (t+1)^num_within_levels);
        outcomes = 1:(t+1)^num_within_levels;
        
        %calculate measurement variability for each subject based upon
        %condition means and number of trials
        for n = 1:nPilotSubs(between_level)
            
            cond_probs = cell(1, num_within_levels);
            
            for c = 1:num_within_levels
                cond_probs{between_level, c} = binopdf(0:t, t, prefs.Y(prefs.subs == between_subs{between_level}(n) & prefs.F1 == between_level & prefs.F2 == c));
            end
            
            %first two conditions
            tmp =  cond_probs{between_level, 1}' * cond_probs{between_level, 2};
            tmp = reshape(tmp, numel(tmp), 1);
            
            %additional conditions
            if num_within_levels > 2
                for c = 3:num_within_levels
                    tmp = reshape(tmp*cond_probs{between_level, c}, numel(tmp)*length(cond_probs{between_level, c}), 1);
                end
            end
            
            %cond_score = cell(1, nConds);
            for c = 1:num_within_levels
                cond_score{between_level, c} = repmat(repelem(0:t, (t+1)^(c-1))', (t+1)^(num_within_levels-c), 1);
            end
            
            cFinal{between_level}(n, :) = tmp;
            
        end
    end
    
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        pc = round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)));
        disp([num2str(pc), '% Complete']);
        
        nSubs_Total = sub_vector(sub_count);
        %nSubs = zeros(1,nConds);
        
        num_trials(trial_count, sub_count) = t;
        requested_sample_size(trial_count, sub_count) = nSubs_Total;
        
        for between_level = 1:num_between_levels
            %number of subjects to simulate
            nSubs(between_level) = round(nSubs_Total*prefs.condition_allocation(between_level));
            outcome_samples{between_level} = randsample(outcomes, nSims*nSubs(between_level), 'true', mean(cFinal{between_level}));
        end
        
        sample_size(trial_count, sub_count) = sum(nSubs);
        
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
                diff_scores = cond_score{bl, wl1}(outcome_samples{bl}) - cond_score{bl, wl2}(outcome_samples{bl});
                diff_scores = reshape(diff_scores, nSubs(bl), nSims);
                dz_vect = mean(diff_scores)./std(diff_scores);
                [~,p] = ttest(diff_scores);
                power_marker(comp, :) = p < prefs.alpha & dz_vect > 0;
                
            else %between comparison
                bl1 = level_between(prefs.comps(comp,1));
                bl2 = level_between(prefs.comps(comp,2));
                wl1 = level_within(prefs.comps(comp,1));
                wl2 = level_within(prefs.comps(comp,2));
                
                c1 = reshape(cond_score{bl1, wl1}(outcome_samples{bl1}), nSubs(bl1), nSims);
                c2 = reshape(cond_score{bl2, wl2}(outcome_samples{bl2}), nSubs(bl2), nSims);
                ds_vect = (mean(c1) - mean(c2)) ./...
                    (((nSubs(bl1) - 1)*(std(c1).^2) + (nSubs(bl2) - 1)*(std(c2).^2))/(nSubs(bl1) + nSubs(bl2) - 2)).^.5;
                [~,p] = ttest2(c1,c2);
                power_marker(comp, :) = p < prefs.alpha & ds_vect > 0;
                
            end
            
            
            
        end
        
        if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %%%
            
            %anova part
            Y = [];
            S = [];
            WF = [];
            BF = [];
            sm = 1;
            
            for between_level = 1:num_between_levels
                
                for c = 1:num_within_levels
                    Y = [Y;reshape(cond_score{between_level, c}(outcome_samples{between_level}), nSubs(between_level), nSims)];
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

end

