function options = get_gp_default_options
options.modelType = 'GPML';
options.pca = 0;
options.logModel = 1;
options.ppSize = 300;
options.trainSubSize = 300;
options.doQuadratic = 0;

options.internal_opt = 1;
options.hyp_opt_obj = 'mll';
options.hyp_opt_algorithm = 'minFunc';
options.hyp_opt_steps = 50;
options.paramsLowerBound = -3;
options.paramsUpperBound = 3;

options.unique_model_name = 'GP';