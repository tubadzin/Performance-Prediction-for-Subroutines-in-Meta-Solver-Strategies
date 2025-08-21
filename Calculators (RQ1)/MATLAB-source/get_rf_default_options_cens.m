function options = get_rf_default_options_cens

options = cens_get_rf_default_options;
options.censoring = 'fill_in_mean';
options.unique_model_name = 'RF-def-fill_in_mean';