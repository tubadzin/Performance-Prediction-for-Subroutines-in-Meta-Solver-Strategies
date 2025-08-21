function options = fill_in_options(options)
options.remove_constant_Theta = 1;
options.ignore_conditionals = 1;
options.kappa_max = 300;
options.overallobj = 'mean';
options.cutoff_penalty_factor = 1;
options.min_variance = 1e-2; % 1e-14; % 1e-2 for modelling paper, better LL
