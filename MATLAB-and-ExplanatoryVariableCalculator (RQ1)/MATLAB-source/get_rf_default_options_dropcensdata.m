function options = get_rf_default_options_dropcensdata

options = cens_get_rf_default_options;
options.censoring = 'dropcensdata';
options.unique_model_name = 'RF-def-dropcensdata';