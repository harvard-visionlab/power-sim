function [f1, f2] = Cond2FactorLevels(cond, prefs)
%FACTORLEVELS2COND convert factor levels to condition number

total_conditions = prefs.f1_num_levels * prefs.f2_num_levels;
f2 = rem(cond, prefs.f2_num_levels);
if f2 == 0
    f2 = prefs.f2_num_levels;
end

if cond == total_conditions
    f1 = prefs.f1_num_levels;
else
    f1 = ceil(rem(cond, total_conditions)/prefs.f2_num_levels);
end
    

end

