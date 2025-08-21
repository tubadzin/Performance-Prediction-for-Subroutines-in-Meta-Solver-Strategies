function options = get_rf_default_options_ignorecensinfo

options = cens_get_rf_default_options;
options.censoring = 'ignorecensinfo';
options.unique_model_name = 'RF-def-ignorecensinfo';