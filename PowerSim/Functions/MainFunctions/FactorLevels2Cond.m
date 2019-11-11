function cond = FactorLevels2Cond(f1,f2, prefs)
%FACTORLEVELS2COND convert factor levels to condition number

cond = prefs.f1_num_levels*(f1-1) + f2;
end

