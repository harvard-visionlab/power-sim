function prefs = MainSetup(prefs)
prefs.max_array_size = 1000000;
prefs.figure_width = 700;
prefs.figure_height = 700;

%randomization
try
    rng('default');
    rng('shuffle');
catch
    matlabVersion = version('-release');
    if str2num(matlabVersion(1:4)) < 2009
        rand  ('twister', sum(100*clock));
    else
        RandStream.setDefaultStream(RandStream('mt19937ar','seed', sum(100*clock)));
    end
end
end
