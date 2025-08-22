function epm_experiments(idx,exps)
% Carries out the experiment for simple empirical hardness models
% (i.e., only features vary) for the distribution with index idx.

if isdeployed
    idx = str2double(idx); % 1 would be just the bigmix with cplex csv
    exps = str2double(exps); % 1 is :Evaluate performance by cross validation.
end

if ismac
    dir = '/Users/piotrmalkowski/Desktop/Bachelor/Hutter et al./data/';
elseif isunix
    dir = 'data/';
else
    dir = 'data\';
end

distributions = {};

% % MIP
distributions{end+1} = {'MIP/BIGMIX-train_test-features-withfilename.csv', 'MIP/BIGMIX-cplex.csv', 1e7};
distributions{end+1} = {'MIP/BIGMIX-train_test-features-withfilename.csv', 'MIP/BIGMIX-gurobi.csv', 1e7};
distributions{end+1} = {'MIP/BIGMIX-train_test-features-withfilename.csv', 'MIP/BIGMIX-scip.csv', 1e7};
distributions{end+1} = {'MIP/BIGMIX-train_test-features-withfilename.csv', 'MIP/BIGMIX-lpsolve.csv', 1e7};

distributions{end+1} = {'MIP/CORLAT-train_test-features-withfilename.csv', 'MIP/CORLAT-cplex.csv'};
distributions{end+1} = {'MIP/CORLAT-train_test-features-withfilename.csv', 'MIP/CORLAT-gurobi.csv'};
distributions{end+1} = {'MIP/CORLAT-train_test-features-withfilename.csv', 'MIP/CORLAT-scip.csv'};
distributions{end+1} = {'MIP/CORLAT-train_test-features-withfilename.csv', 'MIP/CORLAT-lpsolve.csv'};

distributions{end+1} = {'MIP/RCW-train_test-features-withfilename.csv', 'MIP/RCW-cplex.csv'};
distributions{end+1} = {'MIP/REG-train_test-features-withfilename.csv', 'MIP/REG-cplex.csv'};
distributions{end+1} = {'MIP/CORLAT-REG-features.csv', 'MIP/CORLAT-REG-cplex.csv'}; 
distributions{end+1} = {'MIP/CORLAT-REG-RCW-features.csv', 'MIP/CORLAT-REG-RCW-cplex.csv'}; 


distributions{end+1} = {'dummy-feats.csv', 'dummy-runtime.csv'};

do_exp = 1;
do_plots = 1;
num_folds = 10;

% FOR RUNNING EXPERIMENTS:
% All models we'll ever need (except RR-elim).
% options_vec = {get_lr_default_options,get_spore_options_rmse,get_lr_jacm_options,get_nn_options,get_rf_default_options, get_lr_cv_options,get_spore_cv_options_rmse,get_nn_cv_options,get_gp_default_options,get_regtree_default_options,get_rf_cv_options};

% % RR-elim (takes forever).
% % options_vec = {get_lr_jacm_options};
% 
% % FOR PLOTS (set stats to display in plot_ehm_results.m):
% % First bake-off:
% %options_vec = {get_lr_cv_options,get_lr_jacm_options,get_spore_cv_options_rmse,get_nn_cv_options,get_gp_default_options,get_regtree_default_options,get_rf_default_options,get_rf_cv_options};
% %options_vec = {get_lr_cv_options,get_lr_jacm_options,get_spore_cv_options_rmse,get_nn_cv_options,get_gp_default_options,get_regtree_default_options,get_rf_default_options};
% 
% % First bake-off with defaults:
% options_vec = {get_lr_default_options,get_lr_jacm_options,get_spore_options_rmse,get_nn_options,get_gp_default_options,get_regtree_default_options,get_rf_default_options};
% options_vec = {get_lr_default_options,get_lr_jacm_options,get_spore_options_rmse,get_nn_options,get_gp_default_options,get_regtree_default_options,get_rf_default_options};

% % Default vs  hyperparameter optimized version
% % options_vec = {get_lr_default_options,get_lr_cv_options,get_spore_options_rmse, get_spore_cv_options_rmse, get_nn_options, get_nn_cv_options,get_rf_default_options,get_rf_cv_options};%,get_gp_default_options, get_regtree_default_options};
% % 
% % % C/V version, ordering for tables
% % options_vec = {get_lr_cv_options,get_spore_cv_options_rmse,get_nn_cv_options,get_gp_default_options,get_regtree_default_options,get_rf_default_options}%,get_rf_cv_options}

% C/V version, ordering for scaling experiments:
% options_vec = {get_rf_cv_options,get_regtree_default_options,get_gp_default_options,get_nn_cv_options,get_spore_cv_options_rmse,get_lr_cv_options}%,get_rf_cv_options}
% options_vec = {get_regtree_default_options}%,get_rf_cv_options}

% % options_vec = {get_lr_jacm_options}

options_vec = {get_rf_default_options};
if 1 % reverse order for plotting
    tmp = options_vec;
    options_vec = {};
    for i=length(tmp):-1:1
        options_vec{end+1} = tmp{i};
    end
end

cat = [];
catDomains = [];
algo_deterministic = 0;

domain_names = {};
for dist_idx=idx
    feature_filename = strcat(dir, distributions{dist_idx}{1});
    runtime_filename = strcat(dir, distributions{dist_idx}{2});

    domain_name = regexprep(runtime_filename, '\.', '_');
    domain_name = regexprep(domain_name, '\\', '_');
    domain_name = regexprep(domain_name, '/', '_');
    domain_name = regexprep(domain_name, ' ', '_');

    domain_names{end+1} = domain_name;

   if do_exp
        %=== Read instance names and runtimes from runtime file.
        [instance_names] = textread(runtime_filename,'%s%*[^\n]', 'bufsize', 100000);
        for i=1:length(instance_names)
            instance_names{i} = instance_names{i}(1:end-1);
        end
        y = csvread(runtime_filename, 0, 1);

        %=== Read instance names and features from feature file.
        [instance_names_feat_file] = textread(feature_filename,'%s%*[^\n]', 'bufsize', 100000, 'delimiter',',');
        instance_names_feat_file = instance_names_feat_file(2:end);
        all_features = csvread(feature_filename, 1, 1);
        featureNames = textread(feature_filename, '%s', 1, 'whitespace', '\n', 'bufsize', 10000);
        featureNames = strread(featureNames{1},'%s','whitespace',',');
        featureNames = deblank(featureNames);
        featureNames= featureNames(2:end);

        %=== Assert order of instances is the same.
        assert( length(instance_names_feat_file) == length(instance_names) );
        for i=1:length(instance_names_feat_file)
            if ~( strcmp(instance_names_feat_file{i}, instance_names{i}) );
                instance_names_feat_file{i}
                instance_names{i}
                error
            end
        end

sub_idx = 1:size(all_features,1); 
feat_idx = 1:size(all_features,2);

%sub_idx = 1:100;
%feat_idx = [1:30];
all_features = all_features(sub_idx,feat_idx);
y = y(sub_idx);

        axmax = 3600;
        if length(distributions{dist_idx}) > 2
            axmax = distributions{dist_idx}{3};
        end
        for exp = exps
            do_ehm_experiments(options_vec, all_features,y,featureNames,exp,domain_name,cat,catDomains,algo_deterministic,0.001,axmax,num_folds);
        end
   end
end

if do_plots
    for exp = exps
        for dom_num=1:length(domain_names)
            domain = domain_names{dom_num};
            plot_ehm_results(options_vec, domain, exp, dom_num);
        end
    end
end

