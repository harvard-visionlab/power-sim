function power_results = PowerAnalysis_tTest_Within_Fast(prefs)

%set up data (mean scores, nSubs * nConds)
subs = unique(prefs.data(:,1));
nConds = length(unique(prefs.data(:,3)));
nPilotSubs = length(subs); %how many subjects in actual data

for s = 1:nPilotSubs
    for c = 1:nConds
        pilot_data(s, c) = mean(prefs.data(prefs.data(:,1) == subs(s) & ...
            prefs.data(:,3) == c,2));
    end
end

%simulation info
nSims = prefs.nSims; %number of experiments to simulate
sub_vector = prefs.N_range; %number of subs per simulation
trial_vector = fliplr(prefs.trial_range); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

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
            cond_probs{c} = binopdf(0:t, t, pilot_data(n,c));
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
        disp([num2str(round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)))), '% Complete']);
        
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
