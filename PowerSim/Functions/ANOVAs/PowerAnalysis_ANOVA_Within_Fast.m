function power_results = PowerAnalysis_ANOVA_Within_Fast(prefs)

%simulation info
nSims = prefs.nSims; %number of experiments to simulate
nPilotSubs = length(unique(prefs.data(:,1))); %how many subjects in actual data
sub_vector = prefs.N_range; %number of subs per simulation
trial_vector = fliplr(prefs.trial_range); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest
nConds = prefs.f1_num_levels*prefs.f2_num_levels; %number of conditions

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = zeros(length(trial_vector), length(sub_vector));
dz_vect = cell(1, nComps);

for trial_count= 1:length(trial_vector)
    
    t = trial_vector(trial_count);
    
    %determine condition difference pdf for each subject in pilot data
    cFinal = zeros(nPilotSubs, (t+1)^nConds);
    outcomes = 1:(t+1)^nConds;
    
    %calculate measurement variability for each subject based upon
    %condition means and number of trials
    for n = 1:nPilotSubs
        cond_probs = cell(1, nConds);
        for c = 1:nConds
            cond_probs{c} = binopdf(0:t, t, prefs.sub_means(n,c));
        end
        
        %first two conditions
        tmp =  cond_probs{1}' * cond_probs{2};
        tmp = reshape(tmp, numel(tmp), 1);
        
        %additional conditions
        if nConds > 2
            for c = 3:nConds
                tmp = reshape(tmp*cond_probs{c}, numel(tmp)*length(cond_probs{c}), 1);
            end
        end
        
        cond_score = cell(1, nConds);
        
        for c = 1:nConds
            cond_score{c} = repmat(repelem(0:t, (t+1)^(c-1))', (t+1)^(nConds-c), 1);
        end
   
        cFinal(n, :) = tmp;
    end
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        pc = round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)));
        disp([num2str(pc), '% Complete']);
        
        %number of subjects to simulate
        nSubs = sub_vector(sub_count);
        outcome_samples = randsample(outcomes, nSims*nSubs, 'true', mean(cFinal));
        
        sample_size(trial_count, sub_count) = nSubs;
        num_trials(trial_count, sub_count) = t;
        
        % do condition comparisons
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            diff_scores = cond_score{prefs.comps(comp,1)}(outcome_samples) - cond_score{prefs.comps(comp,2)}(outcome_samples);
            diff_scores = reshape(diff_scores, nSubs, nSims);
            dz_vect{comp}{trial_count, sub_count} = mean(diff_scores)./std(diff_scores);
            [~,p] = ttest(diff_scores);
            power_marker(comp, :) = p < prefs.alpha & dz_vect{comp}{trial_count, sub_count} > 0;
        end
        
        if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %anova part
            Y = [];
            for c = 1:nConds
                Y = [Y;reshape(cond_score{c}(outcome_samples), nSubs, nSims)];
            end
            S = repmat(1:nSubs, 1, nConds)';
            
            F1 = [ones(nSubs * nConds/2, 1); 1 + ones(nSubs * nConds/2, 1)];
            F2 = repmat([ones(nSubs, nConds/4); 1 + ones(nSubs, nConds/4)], 2, 1);
            
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


end