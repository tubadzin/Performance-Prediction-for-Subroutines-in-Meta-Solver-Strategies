function options = get_rf_default_options_fill_in_samples

options = cens_get_rf_default_options;
options.censoring = 'fill_in_samples';
options.unique_model_name = 'RF-def-fill_in_samples';