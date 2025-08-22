function do_ehm_experiments(options_vec,X,y,featureNames,exp,domain,cat,catDomains,algo_deterministic,axmin,axmax,num_folds)
seed = 1234;
if nargin < 12
    num_folds = inf;
end
if nargin < 10
    axmin = 0.001;
    axmax = 3600;
end
% Given data X and y, experiment with different models.
outdomain = strcat('results/ehm/', domain, '/');
mkdir(outdomain);

y = max(y,0.005);
log_y = log10(y);
mean_log_y = mean(log_y);
std_log_y = std(log_y);
median_log_y = median(log_y);
q25_log_y = quantile(log_y,0.25);
q75_log_y = quantile(log_y,0.75);
min_log_y = min(log_y);
max_log_y = max(log_y);
csvwrite(strcat(outdomain, '-stats.csv'), [mean_log_y,std_log_y,median_log_y,q25_log_y,q75_log_y,min_log_y,max_log_y]);

fprintf(strcat([' Statistics of log10-transformed y for domain ', domain, ':\n mean=', num2str(mean_log_y), '\n std = ', num2str(std_log_y), '\n median = ', num2str(median_log_y), '\n q_25 = ', num2str(q25_log_y), '\n q_75 = ', num2str(q75_log_y), '\n min = ', num2str(min_log_y), '\n max = ', num2str(max_log_y), '\n']));


fprintf(strcat(['EHM experiments for domain ', domain, '\n']));




%% Experiment 1: Evaluate performance by cross validation.
if exp==1
    fprintf(strcat(['EHM experiment 1 for domain ', domain, '\n']));
    
    for model_idx=1:length(options_vec)
        rand('twister', seed);
        fprintf(strcat([' Model ', num2str(model_idx), '/', num2str(length(options_vec)), ' (', options_vec{model_idx}.unique_model_name, '):\n']));
        [rmses, lls, ccs, cc_ranks, learnTimes, predTimes, y_to_plot, all_y_cross, all_y_cross_var, cens_to_plot] = crossval_model_perf(X, y, num_folds, options_vec{model_idx}, cat, catDomains, algo_deterministic, 1, domain, axmin, axmax);
        filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_CV-result.mat');
        mkdir(strcat(outdomain, options_vec{model_idx}.unique_model_name));
        fprintf(strcat(['Saving CV results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
        save(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'learnTimes', 'predTimes', 'y_to_plot', 'all_y_cross', 'all_y_cross_var', 'cens_to_plot');
    end
end

%% Experiment 2: scaling performance with # data points
if exp==2
    
    fprintf(strcat(['EHM experiment 2 for domain ', domain, '\n']));

    maxNumTrain = length(y)*9/10;
    numTicks = 8;
    numTrains = [1];
    for i=1:numTicks-1
        numTrains(end+1) = numTrains(end) * maxNumTrain^(1/(numTicks-1));
    end
    for i=1:numTicks
        numTrains(i) = ceil(numTrains(i));
    end
       
    for model_idx = 1:length(options_vec)
        fprintf(strcat([' Model ', num2str(model_idx), '/', num2str(length(options_vec)), ' (', options_vec{model_idx}.unique_model_name, '):\n']));
        for i=1:length(numTrains)
            numTrain = numTrains(i);
            fprintf(strcat(['  numTrain = ', num2str(numTrain), '\n']));

            [rmses, lls, ccs, cc_ranks, timesToLearn, timesToPredict] = repeated_simple_model_perf(X, y, numTrain, options_vec{model_idx}, cat, catDomains, algo_deterministic, 10);

            filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_scaleN', num2str(numTrain), '-result.mat');
            fprintf(strcat(['Saving scaling results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
            save(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'timesToLearn', 'timesToPredict');
        end
    end
end



%% Experiment 3: scaling performance with #trees, size of active set, and #features.
if exp==3
    fprintf(strcat(['EHM experiment 3 for domain ', domain, '\n']));

    numTrain = 1000;
    algo_deterministic = 0;

    for model_idx = 1:length(options_vec)
        options = options_vec{model_idx};
        switch options.modelType
            case 'LR'
                numScal = [1,2,4,8,16,32,64,128];
            case 'rf'
                if strcmp(options.unique_model_name, 'RF-cv')
                    numScal = [1,2,4,8,16,32,64,128];
                else
                    numScal = [1,2,4,8,16,32,64,128,256,512,1024];
                end
            case 'GPML'
                numScal = [1,2,4,8,16,32,64,128,256,512,1024];
            case 'regression-tree'
                numScal = [1];
            case 'spore'
                numScal = [1,2,4,8,16,32,64,128];
            case 'nn'
                numScal = [1,2,4,8,16,32,64,128];
        end
        for i=1:length(numScal)
            fprintf(strcat(['  Version ', num2str(i), '/', num2str(length(numScal)), ' ...\n']));
            switch options.modelType
                case 'LR'
                    maxModelSize = numScal(i);
                    maxModelSize = (log(maxModelSize)-log(options.maxModelSize_min)) / (log(options.maxModelSize_max)-log(options.maxModelSize_min));
                    options.tuning_params(1) = maxModelSize;

                case 'rf'
                    options.nSub = numScal(i);
                case 'GPML'
                    options.ppSize = numScal(i);
                    options.trainSubSize = numScal(i);
                case 'regression-tree'
                case 'spore'
                    max_terms = numScal(i);
                    max_terms = (log(max_terms)-log(options.max_terms_min)) / (log(options.max_terms_max)-log(options.max_terms_min)); 
                    options.tuning_params(3) = max_terms;
                    
                case 'nn'
                    nhidden = numScal(i);
                    nhidden = (log(nhidden)-log(options.nhidden_min)) / (log(options.nhidden_max)-log(options.nhidden_min));
                    options.tuning_params(1) = nhidden;
            end           
            [rmses, lls, ccs, cc_ranks, timesToLearn, timesToPredict] = repeated_simple_model_perf(X, y, numTrain, options, cat, catDomains, algo_deterministic, 10);
            
            filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_scaleC', num2str(numScal(i)), '-result.mat');
            fprintf(strcat(['Saving model complexity scaling results for model ', options_vec{model_idx}.unique_model_name, ' to file ', filename, '\n']));
            save(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'timesToLearn', 'timesToPredict');
        end
    end
end

