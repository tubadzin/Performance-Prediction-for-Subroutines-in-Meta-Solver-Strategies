function options = get_regtree_default_options
options.modelType = 'regression-tree';
options.logModel = 1;
options.pca = 0; % do no PCA here, for feature importance
options.init_Splitmin = 5;
options.pruneRatio = 0.5;
options.strategyForMissing = 0;
options.doQuadratic = 0;

options.unique_model_name = 'RT';