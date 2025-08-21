This is the source code for running the experiments for our AIJ submission "Algorithm runtime prediction: The New State of the Art".
The main file for experiments in the features space is epm_experiments.m
The main file for experiments in the configuration space and the joint feature/configruation space is epm_matrix_experiments.m

Both of these files can both run experiments and plot their results, depending on whether
do_exp and/or do_plots is set in the file.

Both of these files take some input parameters to only consider part of the experiments.
In fact, we ran the experiments on a compute cluster, one experiment at a time.
The callstrings we used are in the following files:
featurespace_exps.txt
matrix_exps.txt
censoring_experiment_cmdline_calls.txt

For questions and comments, please contact Frank Hutter, hutter@cs.ubc.ca