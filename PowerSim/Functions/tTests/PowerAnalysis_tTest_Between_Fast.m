function power_results = PowerAnalysis_tTest_Between_Fast(prefs)

if sum(prefs.condition_allocation) ~= 1
    error('Condition allocation does not add up to 1 (100%). Use fractions if this is due to rounding error (e.g., use 1/3 instead of .33)')
elseif length(prefs.condition_allocation) ~= length(prefs.cond_names)
    error('Must have an allocation amount for each between-subjects level.')
end

%Get data setup properly (array slot for each condition, 
%simulation info
nSims = prefs.nSims; %number of experiments to simulate
nConds = length(unique(prefs.data(:,3)));
sub_vector = prefs.N_range; %number of subs per simulation
trial_vector = fliplr(prefs.trial_range); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

for c = 1:nConds
   cond_subs{c} = unique(prefs.data(prefs.data(:,3) == c,1));
   nPilotSubs(c) = length(cond_subs{c}); %how many subjects in actual data per condition
   for s = 1:nPilotSubs(c)
      pilot_data{c}(s) = mean(prefs.data(prefs.data(:,1) == cond_subs{c}(s) & ...
            prefs.data(:,3) == c,2));
   end
end


%preallocate
power = zeros(length(trial_vector), length(sub_vector));
requested_sample_size = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
subs_by_cond = cell(length(trial_vector), length(sub_vector));
num_trials = zeros(length(trial_vector), length(sub_vector));
ds_vect = cell(1, nComps);


for trial_count= 1:length(trial_vector)
    
    t = trial_vector(trial_count);
    
    %determine condition difference pdf for each subject in pilot data
    cFinal = cell(1, nConds);
    outcomes = 0:t;
    
    %calculate measurement variability for each subject based upon 
    %condition means and number of trials
    
    cond_prob = cell(1, nConds);
    for c = 1:nConds
        for n = 1:nPilotSubs(c)
            cond_prob{c}(n,:) = binopdf(0:t, t, pilot_data{c}(n));
        end
        cFinal{c} = mean(cond_prob{c});
    end
    
    
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        disp([num2str(round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)))), '% Complete']);
        
        %number of subjects to simulate per condition
        
        nSubs_Total = sub_vector(sub_count);
        outcome_samples = cell(1, nConds);
        nSubs = zeros(1,nConds);

        for c = 1:nConds
            nSubs(c) = round(nSubs_Total*prefs.condition_allocation(c));
            outcome_samples{c} = randsample(outcomes, nSims*nSubs(c), 'true', cFinal{c});
        end
        
        subs_by_cond{trial_count, sub_count} = nSubs;
        requested_sample_size(trial_count, sub_count) = nSubs_Total;
        sample_size(trial_count, sub_count) = sum(nSubs);
        num_trials(trial_count, sub_count) = t;
        
        % do condition comparisons
        %p = cell(1, nComps);
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            c1 = reshape(outcome_samples{prefs.comps(comp,1)}, nSubs(prefs.comps(comp,1)), nSims);
            c2 = reshape(outcome_samples{prefs.comps(comp,2)}, nSubs(prefs.comps(comp,2)), nSims);
            ds_vect{comp}{trial_count, sub_count} = (mean(c1) - mean(c2)) ./...
                (((nSubs(prefs.comps(comp,1)) - 1)*(std(c1).^2) + (nSubs(prefs.comps(comp,2)) - 1)*(std(c2).^2))/(nSubs(prefs.comps(comp,1)) + nSubs(prefs.comps(comp,2)) - 2)).^.5;
            [~,p] = ttest2(c1,c2);
            power_marker(comp, :) = p < prefs.alpha & ds_vect{comp}{trial_count, sub_count} > 0;
            
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
power_results.ds_vect = ds_vect; %effect size vector for each design
power_results.sub_by_cond = subs_by_cond;  %subs in each cond for each design
power_results.sub_vector = sample_size(1,:);
power_results.trial_vector = trial_vector;
end