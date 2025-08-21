function options = get_rf_default_options
options.modelType = 'rf';
options.nSub = 10;
options.logModel = 1;
options.pca = 0; % do no PCA here, for feature importance
options.orig_rf = 0;
options.strategyForMissing = 0;
options.doQuadratic = 0;

% Tuning parameters:
splitMin_def = 5;
split_ratio_def = 0.5;

options.paramsLowerBound = 0;
options.paramsUpperBound = 1;

options.splitMin_max = 10;
options.splitMin_min = 1;
options.split_ratio_max = 1;
options.split_ratio_min = 0;

Splitmin = (splitMin_def-options.splitMin_min) / (options.splitMin_max-options.splitMin_min);
split_ratio = (split_ratio_def-options.split_ratio_min) / (options.split_ratio_max-options.split_ratio_min);;

options.tuning_params = [split_ratio, Splitmin];
options.unique_model_name = 'RF-def';

options.opt = 0;