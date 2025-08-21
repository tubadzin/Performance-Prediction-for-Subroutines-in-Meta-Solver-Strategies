function model_perf_scenario_clean(options_vec, matrix_filename, filename_with_params_for_matrix, feature_file, tuningScenario, filename_with_config_space, captime_for_matrix, exp, singleExp, capping_type)
% MODEL_PERF_SCENARIO: generates diagnostic plots to study how well
% various regression methods predict performance
%
% Input: tuningScenario: also used as prefix for output plots
%        matrix_filename: name of file with I lines for I instances, each
%        line <instance_name>, <performance of config1>, <performance of
%        config2>, ... In my experiments: config1 = default config (first
%        to be explored in IS-best)
%        captime_for_matrix: time after which runs for generating the
%        matrix were terminated -- matrix entries for terminated runs must
%        have values larger than this captime
%        filename_with_params_for_matrix: name of file with C parameter 
%        configurations, where C is the number of columns in the input 
%        matrix, and the results in column i are for configuration i


isTinyTest = 0;
tuningScenario

%=== Read in a representation of the "function" we optimize, including info
%=== on parameters.
func = read_func(tuningScenario, filename_with_config_space);

%=== Read in matrix data.
data = csvread(matrix_filename,0,1);
matrix = data';
[C,I] = size(matrix); % C: #parameter configurations; I: #instances

%=== Prepare data using PAR-1.
cutoff_multiplier = 1;
for i=1:size(matrix,2)
    matrix(find(matrix(:,i)>captime_for_matrix-0.01), i) = captime_for_matrix * cutoff_multiplier; 
end
matrix = max(matrix, 0.005);

cpu_days_for_data_gathering = sum(sum(matrix))/3600/24

%=== Read in configurations used to generate the matrix data.
Theta = read_param_configs(func, filename_with_params_for_matrix, 1);
assert(all(all(Theta(:,func.cat) >= 1)));

if exp==0
    def_perf = matrix(1,:);
    
    idxs = find(def_perf>10);
    [tmp, idx] = min(def_perf(idxs));
    chosen_idx = idxs(idx);
    if isempty(chosen_idx)
        [tmp, chosen_idx] = max(def_perf);
    end
    
    X = Theta;
    y = matrix(:,chosen_idx);
    featureNames = func.param_names;
    cat = func.cat;
    catDomains = func.all_values;
    
    do_ehm_experiments(options_vec, X,y,featureNames,singleExp,strcat('params-',tuningScenario),cat,catDomains,0);
    return
end

%=== Read instance names from matrix data.
[instance_names] = textread(matrix_filename,'%s%*[^\n]', 'bufsize', 100000, 'delimiter', ',');
for i=1:length(instance_names)
    instance_names{i} = instance_names{i}(1:end-1);
end

%=== Read instance names and features from feature file.
all_features = csvread(feature_file, 1, 1);
xNames = textread(feature_file, '%s', 1, 'whitespace', '\n', 'bufsize', 10000);
xNames = strread(xNames{1},'%s','whitespace',',');
xNames = deblank(xNames);
xNames= xNames(2:end);
[instance_names_feat_file] = textread(feature_file,'%s%*[^\n]', 'bufsize', 100000, 'delimiter', ',');
instance_names_feat_file = instance_names_feat_file(2:end);
for i=1:length(instance_names_feat_file)
    instance_names_feat_file{i} = instance_names_feat_file{i}(1:end-1);
end

%=== Assert order of instances is the same.
assert( length(instance_names_feat_file) == length(instance_names) );
for i=1:length(instance_names_feat_file)
    assert( strcmp(instance_names_feat_file{i}, instance_names{i}) );
end

names = [func.param_names; xNames];

%=== Shuffle training and test configurations.
perm = randperm(size(matrix,1));
matrix = matrix(perm,:);
Theta = Theta(perm,:);

%=== Shuffle training and test instances.
perm = randperm(size(matrix,2));
matrix = matrix(:,perm);
all_features = all_features(perm,:);

outdomain = strcat('results/matrix/', tuningScenario, '/');

%% After this setup, do some experiments with the matrix.
%% Experiment 1: matrix predictions
if exp==1
%    options_vec = {get_rf_cv_options,get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options};
%     options_vec = {get_gp_default_options};
    seed = 1234;
    nTrain = 100;
    nTrain = 1000;
    nTrain = 5000000;
    nTrain = 10000;
    matrixPlots = 1;
    if isTinyTest
        nTrain = 10;
    end
    for model_idx=1:length(options_vec)
        options = options_vec{model_idx};
        [timeToLearn, timeToPredict, results, model] = do_matrix_exp(func, matrix, Theta, all_features, names, options, seed, matrixPlots, strcat(tuningScenario,'-',num2str(nTrain)), nTrain, 10000);

        filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_matrix_preds.mat');
        mkdir(strcat(outdomain, options_vec{model_idx}.unique_model_name));
        fprintf(strcat(['Saving results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
        save(filename, 'results', 'timeToLearn', 'timeToPredict', 'nTrain');
    end

    
%     fprintf(strcat([tuningScenario, ' & ']));
%     two_output(timeToLearn{1}, timeToLearn{2}, 1);
%     fprintf(strcat(['& RMSE']));
% 
%     for quadr=1:4
%         fprintf(' & ');
%         n_output([results{1}{quadr}(1), results{2}{quadr}(1), results{3}{quadr}(1)], 1);
%     end
%     fprintf('\\\\\n');
% 
%     fprintf(strcat([' & & & CC']));
%     for quadr=1:4
%         fprintf(' & ');
%         n_output([results{1}{quadr}(3), results{2}{quadr}(3), results{3}{quadr}(3)], 0);
%     end
%     fprintf('\\\\\n');
% 
%     fprintf(strcat([' & & & LL']));
%     for quadr=1:4
%         fprintf(' & ');
%         n_output([results{1}{quadr}(2), results{2}{quadr}(2), results{3}{quadr}(2)], 0);
%     end
%     fprintf('\\\\\n');
end


%% Experiment 2: scaling with training data size
if exp==2
    numTrains = [10,10*sqrt(10),100,100*sqrt(10),1000,1000*sqrt(10),10000,10000*sqrt(10),100000];
    if isTinyTest
        numTrains = [10,10*sqrt(10)];
    end
    matrixPlots = 0;
    seeds = 1:10;
    for model_idx=1:length(options_vec)
        for i=1:length(numTrains)
            nTrain = floor(numTrains(i))
            for j=1:length(seeds)
                seed = seeds(j)
                options = options_vec{model_idx};
                try
                    [timeToLearn{i,j}, timeToPredict{i,j}, results{i,j}] = do_matrix_exp(func, matrix, Theta, all_features, names, options, seed, matrixPlots, tuningScenario, nTrain, 10000);
                catch ME % out of memory for GPs.
                    timeToLearn{i,j} = -1;
                    timeToPredict{i,j} = -1;
                    results{i,j} = [inf, -inf, -inf];
                end
            end
            filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '-nTrain', num2str(nTrain), '_matrix-pred-scale-result.mat');
            mkdir(strcat(outdomain, options_vec{model_idx}.unique_model_name));
            fprintf(strcat(['Saving results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
            save(filename, 'results', 'timeToLearn', 'timeToPredict', 'numTrains');
        end
    end
end

% %% Experiment 3: Exp 1 with pca
% if exp==3
%     options_vec = {get_rf_cv_options,get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options};
% %     options_vec = {get_rf_default_options,get_rf_cv_options,get_gp_default_options}; %get_gp_default_options};
% %     options_vec = {get_regtree_default_options,get_rf_default_options,get_rf_cv_options,get_gp_default_options}; %get_gp_default_options};
%     seed = 1234;
%     matrixPlots = 1;
%     nTrain = 10000;
%     if isTinyTest
%         nTrain = 10;
%     end
% %     nTrain = 10000;
%     for model_idx=1:length(options_vec)
%         options = options_vec{model_idx};
%         options.pca = 7;
%         [timeToLearn{model_idx}, timeToPredict{model_idx}, results{model_idx}, model] = do_matrix_exp(func, matrix, Theta, all_features, names, options, seed, matrixPlots, strcat(tuningScenario,'-',num2str(nTrain)), nTrain, 10000);
% 
%         filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_matrix-pred-scale-result.mat');
%         mkdir(strcat(outdomain, options_vec{model_idx}.unique_model_name));
%         fprintf(strcat(['Saving results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
%         save(filename, 'results', 'timeToLearn', 'timeToPredict', 'numTrains');
%     end
% 
% 	savefile = strcat(outdomain, '-matrix-pred-pca');
% 	save(savefile, 'results', 'timeToLearn', 'timeToPredict', 'options_vec');
% end

% %% Experiment 4: predictions of a single configuration on new instances
% if exp==4
%     numSubC = 100;
%     cat = [];
%     catDomains = [];
%     options_vec = {get_rf_cv_options,get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options};
% %     options_vec = {get_rf_cv_options,get_rf_default_options,get_gp_default_options,get_lr_default_options};
%     % options_vec = {get_rf_cv_options,get_rf_default_options};
% 	model_names = [];
%     for i=1:length(options_vec)
%         model_names{end+1} = options_vec{i}.modelType;
%     end
%     
%     trainI = 1:ceil(I/2.0);
%     testI = ceil(I/2.0)+1:I;
%     
%     model_cc_c = -ones(length(options_vec),1);
%     model_cc_all = -ones(length(options_vec),1);
%     
%     perm=randperm(C);
%     subC = perm(1:numSubC);
% 
%     for model_idx=1:length(options_vec)
%         fprintf(strcat([' Model ', num2str(model_idx), '/', num2str(length(options_vec)), '...\n']));
%         
%         ccs_c = [];
%         lls_c = [];
%         rmses_c = [];
%         timesToLearn_c = [];
%         ccs_all = [];
%         for c=subC
%             [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, c, c, trainI, testI);
%             ccs_c(end+1) = cc;
%             lls_c(end+1) = ll;
%             rmses_c(end+1) = rmse;
%             timesToLearn_c(end+1) = timeToLearn;
%             
% % %            [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, 1:C, c, trainI, testI);
% % %            ccs_all(end+1) = cc;
% %             perm=randperm(C);
% %             perm = setdiff(perm, c);
% % 
% %             [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, [c,perm(1:9)], c, trainI, testI);
% %             ccs_all(end+1) = cc;
%         end
%         model_cc_c(model_idx) = mean(ccs_c);
%         model_cc_all(model_idx) = mean(ccs_all);
%         ccs(model_idx) = mean(ccs_c);
%         lls(model_idx) = mean(lls_c);
%         rmses(model_idx) = mean(rmses_c);
%         timesToLearn(model_idx) = mean(timesToLearn_c);
%     end
%     
%     csvwrite(strcat(outdomain, 'pred-of-inst-rmses'), rmses);
%     csvwrite(strcat(outdomain, 'pred-of-inst-ccs'), ccs);
%     csvwrite(strcat(outdomain, 'pred-of-inst-lls'), lls);
%     csvwrite(strcat(outdomain, 'pred-of-inst-timesToLearn'), timesToLearn);
%         
% % 	fprintf(strcat(['Scenario & ', model_names{1}, ' & ', model_names{2}, ' & ', model_names{3}, '\\\\\n']));
% % 	fprintf(strcat([tuningScenario, ' & ']));
% %     n_output(model_cc_c,0);
% %     fprintf('\\\\\n');
%     modelNames=cell(length(options_vec),1);
%     for i=1:length(modelNames)
%         modelNames{i,1} = options_vec{i}.modelType;
%     end
%     stats_names = {'RMSE', 'CC', 'LL', 'Time to learn'};
%     header_output(modelNames, stats_names);
%     fprintf(strcat(tuningScenario, ' & '));
%     n_output(rmses,1);
%     n_output(ccs,0);
%     n_output(lls,0);
%     n_output(timesToLearn,1);
%     fprintf('\n');    
%     
%     
% % 	fprintf(strcat(['Scenario & Training subset & ', model_names{1}, ' & ', model_names{2}, ' & ', model_names{3}, '\\\\\n']));
% % 	fprintf(strcat([tuningScenario, ' & Single $\\vtheta$ & ']));
% %     n_output(model_cc_c,0);
% %     fprintf('\\\\\n');
% %     fprintf(strcat([tuningScenario, ' & Multiple $\\vtheta \\in \\vTheta$ & ']));
% %     n_output(model_cc_all,0);
% %     fprintf('\\\\\n');
% end

% %% Experiment 2a: scaling performance with # data points
% if exp5
%     options_vec = {get_rf_default_options,get_gp_default_options};%,get_lr_default_options};
% 
%     numTrains = [2,4,8,16,32,64,128,256,512];
%     for i=1:length(numTrains)
%         numTrain = numTrains(i);
%         [rmses, lls, ccs, timesToLearn, timesToPredict] = get_inst_pred_perf(options_vec, matrix, Theta, all_features, C, I, numTrain)
%         for model_idx=1:length(options_vec)
%             rmse_avg(i,model_idx) = mean(rmses(model_idx,:));
%             rmse_std(i,model_idx) = std(rmses(model_idx,:));
%             lls_avg(i,model_idx) = mean(lls(model_idx,:));
%             lls_std(i,model_idx) = std(lls(model_idx,:));
%             cc_avg(i,model_idx) = mean(ccs(model_idx,:));
%             cc_std(i,model_idx) = std(ccs(model_idx,:));
%             timesToLearn_avg(i,model_idx) = mean(timesToLearn(model_idx,:));
%             timesToLearn_std(i,model_idx) = std(timesToLearn(model_idx,:));
%             timesToPredict_avg(i,model_idx) = mean(timesToPredict(model_idx,:));
%             timesToPredict_std(i,model_idx) = std(timesToPredict(model_idx,:));
%         end
%     end
% 
%     model_names = {'RF', 'PP'};
%     plot_Scaling_N(numTrains, cc_avg, cc_std, options_vec, 'Correlation coefficient', model_names, tuningScenario);
%     plot_Scaling_N(numTrains, rmse_avg, rmse_std, options_vec, 'RMSE', model_names, tuningScenario, 'NorthEast');
%     plot_Scaling_N(numTrains, lls_avg, lls_std, options_vec, 'Test Data Log likelihood', model_names, tuningScenario);
%     plot_Scaling_N(numTrains, timesToLearn_avg, timesToLearn_std, options_vec, 'Learn Time', model_names, tuningScenario);
%     plot_Scaling_N(numTrains, timesToPredict_avg, timesToPredict_std, options_vec, 'Prediction Time', model_names, tuningScenario);
% end

% %% Experiment 5: predictions of a single instance with new configurations
% if exp==5
%     numSubI = 100;
%     cat = func.cat;
%     catDomains = func.all_values;
%        
% options_vec = {get_rf_cv_options,get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options};
% % options_vec = {get_rf_cv_options,get_rf_default_options,get_gp_default_options};
% 	model_names = [];
%     for i=1:length(options_vec)
%         model_names{end+1} = options_vec{i}.modelType;
%     end
%     
%     trainC = 1:ceil(C/2.0);
%     testC = ceil(C/2.0)+1:C;
%     
%     model_cc_i = -ones(length(options_vec),1);
%     model_cc_all = -ones(length(options_vec),1);
%     
%     perm=randperm(I);
%     subI = perm(1:numSubI);
% 
%     for model_idx=1:length(options_vec)
%         fprintf(strcat([' Model ', num2str(model_idx), '/', num2str(length(options_vec)), '...\n']));
%         
%         ccs_i = [];
%         lls_i = [];
%         rmses_i = [];
%         timesToLearn_i = [];
%         
%         ccs_all = [];
%         for i=subI
%             [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, trainC, testC, i, i);
%             ccs_i(end+1) = cc;
%             lls_i(end+1) = ll;
%             rmses_i(end+1) = rmse;
%             timesToLearn_i(end+1) = timeToLearn;            
% % %            [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, 1:C, c, trainI, testI);
% % %            ccs_all(end+1) = cc;
% %             perm=randperm(C);
% %             perm = setdiff(perm, c);
% % 
% %             [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options_vec{model_idx}, cat, catDomains, matrix, Theta, all_features, [c,perm(1:9)], c, trainI, testI);
% %             ccs_all(end+1) = cc;
%         end
%         model_cc_i(model_idx) = mean(ccs_i);
%         model_cc_all(model_idx) = mean(ccs_all);
%         ccs(model_idx) = mean(ccs_i);
%         lls(model_idx) = mean(lls_i);
%         rmses(model_idx) = mean(rmses_i);
%         timesToLearn(model_idx) = mean(timesToLearn_i);
%     end
%     
% % 	fprintf(strcat(['Scenario & ', model_names{1}, ' & ', model_names{2}, '\\\\\n']));
% % 	fprintf(strcat([tuningScenario, ' & ']));
% %     n_output(model_cc_i,0);
% %     fprintf('\\\\\n');
% 
%     csvwrite(strcat(outdomain, 'pred-of-conf-rmses'), rmses);
%     csvwrite(strcat(outdomain, 'pred-of-conf-ccs'), ccs);
%     csvwrite(strcat(outdomain, 'pred-of-conf-lls'), lls);
%     csvwrite(strcat(outdomain, 'pred-of-conf-timesToLearn'), timesToLearn);
% 
%     modelNames=cell(length(options_vec),1);
%     for i=1:length(modelNames)
%         modelNames{i,1} = options_vec{i}.modelType;
%     end
%     stats_names = {'RMSE', 'CC', 'LL', 'Time to learn'};
%     header_output(modelNames, stats_names);
%     fprintf(strcat(tuningScenario, ' & '));
%     n_output(rmses,1);
%     n_output(ccs,0);
%     n_output(lls,0);
%     n_output(timesToLearn,1);
%     fprintf('\n');  
% 
% % 	fprintf(strcat(['Scenario & Training subset & ', model_names{1}, ' & ', model_names{2}, ' & ', model_names{3}, '\\\\\n']));
% % 	fprintf(strcat([tuningScenario, ' & Single $\\vtheta$ & ']));
% %     n_output(model_cc_c,0);
% %     fprintf('\\\\\n');
% %     fprintf(strcat([tuningScenario, ' & Multiple $\\vtheta \\in \\vTheta$ & ']));
% %     n_output(model_cc_all,0);
% %     fprintf('\\\\\n');
% end

% if exp==6
%     cat = func.cat;
%     catDomains = func.all_values;
%     [theta_idxs, inst_idxs] = get_idx_pairs(1:C, 1:I);
%     perm = randperm(C*I);
%     nData = 10000;
%     idx = perm(1:nData);
%     theta_idxs = theta_idxs(idx);
%     inst_idxs = inst_idxs(idx);
%     [Xtrain,ytrain] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs);
%     options_vec = {get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options}; 
% %     options_vec = {get_rf_default_options,get_gp_default_options};   
% %     options_vec = {get_regtree_default_options}; %, get_rf_default_options,get_gp_default_options};   
% 
% 
%     if isTinyTest
%         %Tiny data for debugging.
%         p = randperm(size(Xtrain,2));
%         [Xtrain,Xtrain,cat,catDomains,names] = use_feat_subset(Xtrain,Xtrain,cat,catDomains,names,sort(p(1:11)));
%         Xtrain = Xtrain(1:20,:);
%         ytrain = ytrain(1:20);
%     end
%  
%     savefilename = strcat('results/matrix/vimp-sel-val-', outdomain);
%     do_var_imp(options_vec, savefilename, Xtrain, ytrain, cat, catDomains, names);
% end

if exp == 7
%    options_vec = {get_rf_default_options_fill_in_samples,get_rf_default_options_cens,get_rf_default_options_dropcensdata,get_rf_default_options_ignorecensinfo};
%     options_vec = {get_rf_cv_options,get_rf_default_options,get_regtree_default_options,get_gp_default_options,get_lr_default_options};
    switch capping_type
        case 'fixed'
            caps = 10.^((-1:(1+2.4770)/5:2.4770)); % for fixed
        case 'capslack'
            caps = 10.^((0:5)/5); % for capslack
        otherwise
            error
    end
    seeds = 1:10;
    matrixPlots = 1;
    
%     matrixPlots = 1;
    caps = 1;
    seeds = 1:1;
    
    for i=1:length(caps)
        numTrains(i) = 1000;
        nTrain = floor(numTrains(i));
        cap = caps(i);
        for j=1:length(seeds)
            for model_idx=1:length(options_vec)
                seed = seeds(j)
                options = options_vec{model_idx};
%                 try
                   [timeToLearn{model_idx,i,j}, timeToPredict{model_idx,i,j}, results{model_idx,i,j}] = do_matrix_exp(func, matrix, Theta, all_features, names, options, seed, matrixPlots, tuningScenario, nTrain, 10000, cap, capping_type);
%                 catch ME % out of memory for GPs.
%                     timeToLearn{model_idx,i,j} = -1;
%                     timeToPredict{model_idx,i,j} = -1;
%                     results{model_idx,i,j} = [inf, -inf, -inf];
%                 end
%     close all;

            end
        end
    end

	savefile = strcat(outdomain, strcat('cens-matrix-pred-scale-', capping_type));
    mkdir(outdomain);
 	save(savefile, 'results', 'timeToLearn', 'timeToPredict', 'numTrains', 'caps', 'options_vec');
end


function [rmses, lls, ccs, timesToLearn, timesToPredict] = get_inst_pred_perf(options_vec, matrix, Theta, all_features)
[C,I] = size(matrix);
numCruns = min(C,10);
trainI = ceil(I/2.0);
timesToLearn = -ones(length(options_vec), numCruns);
timesToPredict = -ones(length(options_vec), numCruns);
rmses = -ones(length(options_vec), numCruns);
lls = -ones(length(options_vec), numCruns);
ccs = -ones(length(options_vec), numCruns);
perm=randperm(C);
for model_idx=1:length(options_vec)
    fprintf(strcat([' Model ', num2str(model_idx), '/', num2str(length(options_vec)), '...\n']));
    for i=1:numCruns
        fprintf(strcat(['  Configuration ', num2str(i), '/', num2str(numCruns), '...\n']));
        c=perm(i);

        %=== Get indices for training and test instances.
        [theta_idxs, inst_idxs] = get_idx_pairs(c, 1:trainI);
        [Xtrain,ytrain] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs);
        
        [theta_idxs, inst_idxs] = get_idx_pairs(c, trainI+1:I);
        [Xtest,ytest] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs);
        
        % no categorical inputs
        cat = []; 
        catDomains = [];  
        [rmses(model_idx,i), lls(model_idx,i), ccs(model_idx,i), timesToLearn(model_idx,i), timesToPredict(model_idx,i)] = simple_model_perf_train_test(Xtrain, ytrain, Xtest, ytest, cat, catDomains, options_vec{model_idx});

        %=== Call model_perf with right X and y.
%         [rmses(model_idx,i), lls(model_idx,i), ccs(model_idx,i), timesToLearn(model_idx,i), timesToPredict(model_idx,i)] = simple_model_perf(X, y, numTrain, options_vec{model_idx});
    end
end

function [rmse, ll, cc, timeToLearn, timeToPredict] = single_pred_perf(options, cat, catDomains, matrix, Theta, all_features, trainC, testC, trainI, testI)
%=== Get indices for training and test instances.
[theta_idxs, inst_idxs] = get_idx_pairs(trainC, trainI);
maxTrain = 10000;
if length(theta_idxs) > maxTrain
    perm = randperm(length(theta_idxs));
    idx = perm(1:maxTrain);
    theta_idxs = theta_idxs(idx);
    inst_idxs = inst_idxs(idx);
end
[Xtrain,ytrain] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs);

[theta_idxs, inst_idxs] = get_idx_pairs(testC, testI);
[Xtest,ytest] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs);

[rmse, ll, cc, timeToLearn, timeToPredict] = simple_model_perf_train_test(Xtrain, ytrain, Xtest, ytest, cat, catDomains, options);


function func = read_func(tuningScenario, filename_with_config_space)
func.params_filename = filename_with_config_space;
[func.cat, func.cont, func.param_names, func.all_values, func.param_bounds, func.param_trafo, func.is_integer_param, func.default_values, func.cond_params_idxs, func.parent_param_idxs, func.ok_parent_value_idxs] = read_params(func.params_filename);
func.dim = length(func.cat)+length(func.cont);
func.param_names = func.param_names';
func.orig_param_lower_bound = func.param_bounds(:,1);
func.orig_param_upper_bound = func.param_bounds(:,2);
func.param_bounds(find(func.param_bounds(:,1)),1) = 0;
func.param_bounds(find(func.param_bounds(:,2)),2) = 1;
func.default_values = config_transform(func.default_values', func)';

%=== Get number of values for each parameter.
func.num_values = zeros(1,func.dim);
for i=1:func.dim
    func.num_values(i) = length(func.all_values{i});
end
func.deterministic = 0; % necessary for backwards compatibility, but ultimately ignored.