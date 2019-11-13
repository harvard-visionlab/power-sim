function power_results = PowerAnalysis_ANOVA_Between_Fast(prefs)

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


for f1 = 1:prefs.f1_num_levels
    for f2 = 1:prefs.f2_num_levels
        sub_nums{f1,f2} = unique(prefs.sub_nums(prefs.data(:,3) == f1 & prefs.data(:,4) == f2))';
        nPilotSubs(f1, f2) = length(sub_nums{f1,f2});
        for s = 1:nPilotSubs(f1,f2)
            pilot_data{f1,f2}(s) = mean(prefs.data(prefs.data(:,1) == sub_nums{f1,f2}(s) ,2));
        end
    end
end


sub_vector = prefs.N_range; %number of subs per simulation
trial_vector = fliplr(prefs.trial_range); %number of trials per condition
nComps = size(prefs.comps, 1); %number of comparisons of interest

%preallocate
power = zeros(length(trial_vector), length(sub_vector));
requested_sample_size = zeros(length(trial_vector), length(sub_vector));
sample_size = zeros(length(trial_vector), length(sub_vector));
num_trials = zeros(length(trial_vector), length(sub_vector));

for trial_count= 1:length(trial_vector)
    
    t = trial_vector(trial_count);
    
    %determine condition difference pdf for each subject in pilot data
    cFinal = cell(prefs.f1_num_levels, prefs.f2_num_levels);
    outcomes = 0:t;
    
    %calculate measurement variability for each subject based upon 
    %condition means and number of trials
    
    cond_prob = cell(f1, f2);
    for f1 = 1:prefs.f1_num_levels
        for f2 = 1:prefs.f2_num_levels
            for n = 1:nPilotSubs(f1,f2)
                cond_prob{f1,f2}(n,:) = binopdf(0:t, t, pilot_data{f1,f2}(n));
            end
            cFinal{f1,f2} = mean(cond_prob{f1,f2});
        end
    end
     
    %figure out sampling variability based upon number of subjects
    for sub_count = 1:length(sub_vector)
        
        clc
        pc = round(100*((trial_count-1)*length(sub_vector) + sub_count - 1)...
            /(length(trial_vector)*length(sub_vector)));
        disp([num2str(pc), '% Complete']);
        
        nSubs_Total = sub_vector(sub_count);
        outcome_samples = cell(prefs.f1_num_levels, prefs.f2_num_levels);
        nSubs = zeros(prefs.f1_num_levels, prefs.f2_num_levels);
        
        for f1 = 1:prefs.f1_num_levels
            for f2 = 1:prefs.f2_num_levels
                %number of subjects to simulate
                nSubs(f1, f2) = round(nSubs_Total*prefs.condition_allocation(FactorLevels2Cond(f1,f2, prefs)));
                outcome_samples{f1,f2} = randsample(outcomes, nSims*nSubs(f1,f2), 'true', cFinal{f1,f2});
            end
        end
        
        subs_by_cond{trial_count, sub_count} = nSubs;
        requested_sample_size(trial_count, sub_count) = nSubs_Total;
        num_trials(trial_count, sub_count) = t;
        sample_size(trial_count, sub_count) = sum(sum(nSubs));
        
        
        % do condition comparisons
        power_marker = zeros(nComps, nSims);
        for comp = 1:nComps
            
            [c1f1, c1f2] = Cond2FactorLevels(prefs.comps(comp,1), prefs);
            [c2f1, c2f2] = Cond2FactorLevels(prefs.comps(comp,2), prefs);
            
            c1 = reshape(outcome_samples{c1f1, c1f2}, nSubs(c1f1, c1f2), nSims);
            c2 = reshape(outcome_samples{c2f1, c2f2}, nSubs(c2f1, c2f2), nSims);
            ds_vect{comp}{trial_count, sub_count} = (mean(c1) - mean(c2)) ./...
                (((nSubs(c1f1, c1f2) - 1)*(std(c1).^2) + (nSubs(c2f1, c2f2) - 1)*(std(c2).^2))/(nSubs(c1f1, c1f2) + nSubs(c2f1, c2f2) - 2)).^.5;
            [~,p] = ttest2(c1,c2);
            power_marker(comp, :) = p < prefs.alpha & ds_vect{comp}{trial_count, sub_count} > 0;
            
        end
        
        if prefs.sig_ME1 || prefs.sig_ME2 || prefs.sig_int
            %anova part
            Y = [];
            F1 = [];
            F2 = [];
            
            for f1 = 1:prefs.f1_num_levels
                for f2 = 1:prefs.f2_num_levels
                    Y = [Y; reshape(outcome_samples{f1, f2}, nSubs(f1, f2), nSims)];
                    F1 = [F1; repmat(f1, nSubs(f1, f2), 1)];
                    F2 = [F2; repmat(f2, nSubs(f1, f2), 1)];
                end
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
power_results.subs_by_cond = subs_by_cond;  %subs in each cond for each design
power_results.sub_vector = sample_size(1,:);
power_results.trial_vector = trial_vector;

end

