function plot_ehm_results(options_vec, domain, exp, dom_num, axmin, axmax)
if nargin < 6
    axmin = 0.001;
    axmax = 3600+1;
end

outdomain = strcat('results/ehm/', domain, '/');
if exp==1
    %% Experiment 1.
    model_names = {};
    for model_idx = 1:length(options_vec)
        options = options_vec{model_idx};
        model_names{model_idx} = options.unique_model_name;

        filename = strcat(outdomain, options.unique_model_name, '_CV-result.mat');
        load(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'learnTimes', 'predTimes', 'y_to_plot', 'all_y_cross', 'all_y_cross_var', 'cens_to_plot');
        all_rmses(model_idx,:) = rmses;
        all_lls(model_idx,:) = lls;
        all_ccs(model_idx,:) = ccs;
        all_ccranks(model_idx,:) = cc_ranks;
        all_learnTimes(model_idx,:) = learnTimes;
        all_predTimes(model_idx,:) = predTimes;

        figure_prefix = strcat(['results/ehm/', domain, '/', options.unique_model_name]);
        title_prefix = '';
%         plot_simple_pred_scatter(y_to_plot, all_y_cross, all_y_cross_var, cens_to_plot, mean(rmses), mean(ccs), mean(lls), figure_prefix, title_prefix, options.logModel, axmin, axmax, 30);
        close all
    end
%     all_lls = all_lls([4,6,7],:); % only defined for PP and RF
%      all_lls = all_lls([5,7],:); % only defined for PP and RF

% %== For initial comparison
%     stats_names = {'RMSE','Rank correlation coefficient'};
%     if dom_num == 1
%         header_output(model_names, stats_names);
%     end
%     fprintf(strcat(fix_name(fix_domain_name(domain)), ' & '));
%     n_output(mean(all_rmses(1:2,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(3:4,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_ccranks(1:2,:),2),0);   
%     fprintf(' & ');
%     n_output(mean(all_ccranks(3:4,:),2),0);   
%     fprintf('\\\\\n');
%     if ismember(dom_num, [7,11,15,17,21,25,29,32])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end
%     if ismember(dom_num, [17,29])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end      

% %== For hyperparameter optimization
% %     stats_names = {'RMSE','Time to learn [s]'};
%     stats_names = {'Spearman rank correlation coefficient','Log likelihood'};
%     if dom_num == 1
%         header_output(model_names, stats_names);
%     end
%     fprintf(strcat(fix_name(fix_domain_name(domain)), ' & '));
%     n_output(mean(all_ccranks(1:2,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccranks(3:4,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccranks(5:6,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccranks(7:8,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_lls(7:8,:),2),0);
%     
    
%     n_output(mean(all_rmses(1:2,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(3:4,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(5:6,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(7:8,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_learnTimes(1:2,:),2),1);   
%     fprintf(' & ');
%     n_output(mean(all_learnTimes(3:4,:),2),1);   
%     fprintf(' & ');
%     n_output(mean(all_learnTimes(5:6,:),2),1);   
%     fprintf(' & ');
%     n_output(mean(all_learnTimes(7:8,:),2),1);   
%     fprintf('\\\\\n');
%     if ismember(dom_num, [7,11,15,17,21,25,29,32])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end
%     if ismember(dom_num, [17,29])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end      


% additional results for online appendix
%     n_output(mean(all_rmses(1:2,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(3:4,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(5:6,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_rmses(7:8,:),2),1);
%     fprintf(' & ');
%     n_output(mean(all_ccs(1:2,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccs(3:4,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccs(5:6,:),2),0);
%     fprintf(' & ');
%     n_output(mean(all_ccs(7:8,:),2),0);
%     fprintf('\\\\\n');

%     if ismember(dom_num, [7,11,15,17,21,25,29,32])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end
%     if ismember(dom_num, [17,29])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end      

%== RMSE & time.
    stats_names = {'RMSE', 'Time to learn model [s]'};
    if dom_num == 1
        header_output(model_names, stats_names);
    end
    fprintf(strcat(fix_name(fix_domain_name(domain)), ' & '));
    n_output(mean(all_rmses,2),1);
%     n_output(mean(all_ccs,2),0);
    fprintf(' & ');
    n_output(mean(all_learnTimes,2),1);
    fprintf('\\\\\n');
    if ismember(dom_num, [7,11,15,17,21,25,29,32])
        fprintf('\\addlinespace[\\interrowspace]\n');
    end
    if ismember(dom_num, [17,29])
        fprintf('\\addlinespace[\\interrowspace]\n');
    end                    

% %== CC & LL.
%     stats_names = {'Spearman rank correlation coefficient', 'Log likelihood'};
%     if dom_num == 1
%         header_output(model_names, stats_names);
%     end
%     fprintf(strcat(fix_name(fix_domain_name(domain)), ' & '));
%     n_output(mean(all_ccranks,2),0);
%     fprintf(' & ');
%     n_output(mean(all_lls,2),0);
%     fprintf('\\\\\n');
%     if ismember(dom_num, [7,11,15,17,21,25,29,32])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end
%     if ismember(dom_num, [17,29])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end          
    

% %== For remaining evaluations
% %                 stats_names = {'RMSE', 'CC', 'LL', 'Time to learn'};
%     stats_names = {'RMSE', 'Rank correlation coefficient'};
% %     stats_names = {'RMSE', 'Model learning time [s]'};
%     if dom_num == 1
%         header_output(model_names, stats_names);
%     end
%     fprintf(strcat(fix_name(fix_domain_name(domain)), ' & '));
%     n_output(mean(all_rmses,2),1);
% %     n_output(mean(all_ccs,2),0);
%     fprintf(' & ');
% %     n_output(mean(all_lls,2),0);
%     n_output(mean(all_ccranks,2),0);  
% %     n_output(mean(all_learnTimes,2),1);
%     fprintf('\\\\\n');
%     if ismember(dom_num, [7,11,15,17,21,25,29,32])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end
%     if ismember(dom_num, [17,29])
%         fprintf('\\addlinespace[\\interrowspace]\n');
%     end                    
end    

if exp==2
    %% Experiment 2.
    plot_domain = strcat('results/ehm/scale_N_', domain);

    %== Read in results for CV, just to get length of y :-(
    filename = strcat(outdomain, options_vec{1}.unique_model_name, '_CV-result.mat'); 
    load(filename, 'y_to_plot');
    maxNumTrain = length(y_to_plot)*9/10;

    numTicks = 8;
    numTrains = [1];
    for i=1:numTicks-1
        numTrains(end+1) = numTrains(end) * maxNumTrain^(1/(numTicks-1));
    end
    for i=1:numTicks
        numTrains(i) = ceil(numTrains(i));
    end
    numTrains = numTrains(3:end); % don't quite start at 1
    
    model_names = {};
    for model_idx = 1:length(options_vec)
        options = options_vec{model_idx};

        model_names{model_idx} = options.unique_model_name;
        for i=1:length(numTrains)
            numTrain = numTrains(i);
            filename = strcat(outdomain, options.unique_model_name, '_scaleN', num2str(numTrain), '-result.mat');
            try
                load(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'timesToLearn', 'timesToPredict');
                rmse_avg(i,model_idx) = mean(rmses);
                rmse_std(i,model_idx) = std(rmses);
                lls_avg(i,model_idx) = mean(lls);
                lls_std(i,model_idx) = std(lls);
                cc_avg(i,model_idx) = mean(ccs);
                cc_std(i,model_idx) = std(ccs);
                timesToLearn_avg(i,model_idx) = mean(timesToLearn);
                timesToLearn_std(i,model_idx) = std(timesToLearn);
                timesToPredict_avg(i,model_idx) = mean(timesToPredict);
                timesToPredict_std(i,model_idx) = std(timesToPredict);
            catch ME
                rmse_avg(i,model_idx) = inf;
                rmse_std(i,model_idx) = inf;
                lls_avg(i,model_idx) = inf;
                lls_std(i,model_idx) = inf;
                cc_avg(i,model_idx) = inf;
                cc_std(i,model_idx) = inf;
                timesToLearn_avg(i,model_idx) = inf;
                timesToLearn_std(i,model_idx) = inf;
                timesToPredict_avg(i,model_idx) = inf;
                timesToPredict_std(i,model_idx) = inf;     
            end
        end            
    end
    plot_Scaling_N(numTrains, cc_avg, cc_std, options_vec, 'Correlation coefficient', model_names, plot_domain, 'SouthEast');
%     plot_Scaling_N(numTrains, cc_avg, cc_std, options_vec, 'Correlation coefficient', model_names, plot_domain, 'NorthWest');
%     plot_Scaling_N(numTrains, rmse_avg, rmse_std, options_vec, 'RMSE', modelNames, plot_domain, 'NorthEast');
%     plot_Scaling_N(numTrains, lls_avg, lls_std, options_vec, 'Test Data Log likelihood', modelNames, plot_domain);
%     plot_Scaling_N(numTrains, timesToLearn_avg, timesToLearn_std, options_vec, 'Learn Time', modelNames, plot_domain, 'NorthWest', 1);

% %     stats_names = {'RMSE, N=64', 'RMSE, N=512', 'Time, N=64', 'Time, N=512'};
% %     header_output(modelNames, stats_names);
% %     fprintf(strcat(domain, ' & '));
% %     n_output(rmse_avg(5,:),1);
% %     n_output(rmse_avg(8,:),1);
% %     n_output(timesToLearn_avg(5,:),1);
% %     n_output(timesToLearn_avg(8,:),1);
% %     fprintf('\n');    
end

if exp==3
    %% Experiment 3.
    plot_domain = strcat('results/ehm/scale_Compl_by_time_', domain);

    %=== Load results for scaling up model complexity.
    model_names = {};
    namesNumScal = {};
    for model_idx = 1:length(options_vec)
        options = options_vec{model_idx};
        name = fix_modelname(options_vec{model_idx}.unique_model_name);

        switch name
            case 'GP'
                name = 'PP';
            case 'RR'
                name = 'RR(2F)';
            case 'RF-def'
                name = 'RF';
            case 'RT'
                name = 'RT';
            case 'Foba-rmse'
                name = 'RR(FB)';
        end        
        model_names{model_idx} = name;        

        switch options.modelType
            case 'LR'
                numScal = [1,2,4,8,16,32,64,128];
                namesNumScal{model_idx} = '#features';
            case 'rf'
                if strcmp(options.unique_model_name, 'RF-cv')
                    numScal = [1,2,4,8,16,32,64,128];
                else
                    numScal = [1,2,4,8,16,32,64,128,256,512,1024];
                end
                namesNumScal{model_idx} = '#trees';
            case 'GPML'
                numScal = [1,2,4,8,16,32,64];%,128,256,512,1024];
                namesNumScal{model_idx} = 'Size of active set';
            case 'regression-tree'
                numScal = [1];
                namesNumScal{model_idx} = 'dummy';
            case 'spore'
                numScal = [1,2,4,8,16,32,64,128];
                namesNumScal{model_idx} = '#features';
            case 'nn'
                numScal = [1,2,4,8,16,32,64,128];
                namesNumScal{model_idx} = '#neurons';
        end

        for i=1:length(numScal)
            filename = strcat(outdomain, options_vec{model_idx}.unique_model_name, '_scaleC', num2str(numScal(i)), '-result.mat');
            load(filename, 'rmses', 'lls', 'ccs', 'cc_ranks', 'timesToLearn', 'timesToPredict');
            cc_avg{model_idx}(i) = mean(ccs);
            cc_std{model_idx}(i) = std(ccs);
            modelTimesAvg{model_idx}(i) = mean(timesToLearn);
        end
    end
    plot_Scaling_Model_Compl(numScal, namesNumScal, cc_avg, cc_std, modelTimesAvg, 'Correlation coefficient', model_names, plot_domain, 'SouthEast');
end
