function output_tables(domain_names, useRR)
if nargin < 2
    useRR = 1;
end
stats_names = {'RMSE', 'CC', 'LL', 'Time to learn'};
%     options_vec = {get_rf_cv_options,get_rf_default_options,get_gp_default_options,get_lr_default_options};
%     modelNames=cell(length(options_vec),1);

if useRR
    modelNames = {'RF cv', 'RF', 'RT', 'GP', 'RR'};
    indices = [5,4,3,2,1]
else
    modelNames = {'RF cv', 'RF', 'RT', 'GP'};
    indices = [3,4,2,1]
end
modelNames = modelNames(indices);
header_output(modelNames, stats_names);

for i=1:length(domain_names)
    domain = domain_names{i};
    outdomain = strcat('results/ehm/', domain);
    domain = fix_domain_name(domain);

    rmses = csvread(strcat(outdomain, '-all_models-rmses.csv'));
    ccs = csvread(strcat(outdomain, '-all_models-ccs.csv'));
    lls = csvread(strcat(outdomain, '-all_models-lls.csv'));
    learnTimes = csvread(strcat(outdomain, '-all_models-learnTimes.csv'));

    rmses = rmses(indices,:);
    ccs = ccs(indices,:);
    lls = lls(indices,:);
    learnTimes = learnTimes(indices,:);

    fprintf(strcat(fix_name(domain), ' & '));
    n_output(mean(rmses,2),1);
    fprintf(' & ');
    n_output(mean(ccs,2),0);
    fprintf(' & ');
    n_output(mean(lls,2),0);
    fprintf(' & ');
    n_output(mean(learnTimes,2),1);

    fprintf('\\\\\n');    
end

fprintf('\n');    
    
for i=1:length(domain_names)
    domain = domain_names{i};
    outdomain = strcat('results/ehm/', domain);
    domain = fix_domain_name(domain);

    rmses = csvread(strcat(outdomain, '-all_models-rmses.csv'));
    learnTimes = csvread(strcat(outdomain, '-all_models-learnTimes.csv'));

    rmses = rmses(indices,:);
    learnTimes = learnTimes(indices,:);

    fprintf(strcat(fix_name(domain), ' & '));
    n_output(mean(rmses,2),1);
    fprintf(' & ');
    n_output(mean(learnTimes,2),1);

    fprintf('\\\\\n');    
end

fprintf('\n');

for i=1:length(domain_names)
    domain = domain_names{i};
    outdomain = strcat('results/ehm/', domain);
    domain = fix_domain_name(domain);

    ccs = csvread(strcat(outdomain, '-all_models-ccs.csv'));
    lls = csvread(strcat(outdomain, '-all_models-lls.csv'));

    ccs = ccs(indices,:);
    lls = lls(indices,:);

    fprintf(strcat(fix_name(domain), ' & '));
    n_output(mean(ccs,2),0);
    fprintf(' & ');
    n_output(mean(lls,2),0);

    fprintf('\\\\\n');    
end



% fprintf('\n');    
% fprintf('\n');
% %    stats_names = {'CC, N=64', 'CC, N=512', 'Time, N=64', 'Time, N=512'};
% % stats_names = {'CC, N=64', 'CC, N=512'};
% stats_names = {'CC, small N', 'CC, large N'}; % 90\% and less of training data
% header_output(modelNames, stats_names);
% 
% for i=1:length(domain_names)
%     domain = domain_names{i};
%     outdomain = strcat('results/ehm/scale_N/', domain);
%     domain = fix_domain_name(domain);
% 
%     rmse_avg = csvread(strcat(outdomain, '-scale_N_all_models-rmse_avg.csv'));
%     rmse_std = csvread(strcat(outdomain, '-scalN-all_models-rmse_std.csv'));
% 
%     cc_avg = csvread(strcat(outdomain, '-scalN-all_models-cc_avg.csv'));
%     cc_std = csvread(strcat(outdomain, '-scalN-all_models-cc_std.csv'));
% 
%     lls_avg = csvread(strcat(outdomain, '-scalN-all_models-lls_avg.csv'));
%     lls_std = csvread(strcat(outdomain, '-scalN-all_models-lls_std.csv'));
% 
%     timesToLearn_avg = csvread(strcat(outdomain, '-scalN-all_models-timesToLearn_avg.csv'));
%     timesToLearn_std = csvread(strcat(outdomain, '-scalN-all_models-timesToLearn_std.csv'));
% 
%     timesToPredict_avg = csvread(strcat(outdomain, '-scalN-all_models-timesToPredict_avg.csv'));
%     timesToPredict_std = csvread(strcat(outdomain, '-scalN-all_models-timesToPredict_std.csv'));
% 
%     rmse_avg = rmse_avg(:,indices);
%     rmse_std = rmse_std(:,indices);
%     cc_avg = cc_avg(:,indices);
%     cc_std = cc_std(:,indices);
%     lls_avg = lls_avg(:,indices);
%     lls_std = lls_std(:,indices);
%     timesToLearn_avg = timesToLearn_avg(:,indices);
%     timesToLearn_std = timesToLearn_std(:,indices);
%     timesToPredict_avg = timesToPredict_avg(:,indices);
%     timesToPredict_std = timesToPredict_std(:,indices);
% 
%     fprintf(strcat(domain, ' & '));
%     n_output(cc_avg(5,:),0);
%     fprintf(' & ');
%     n_output(cc_avg(8,:),0);
% %        n_output(timesToLearn_avg(5,:),1);
% %        n_output(timesToLearn_avg(8,:),1);
%     fprintf('\\\\\n');    
% end

%     for i=1:length(domain_names)
%         domain = domain_names{i};
%         outdomain = strcat('results/ehm/', domain);
%         
%         for j=1:5
%             rmse_avg{j} = csvread(strcat(outdomain, '-scalComp-model,', num2str(model_idx), '-rmse_avg.csv'));
%             rmse_std{j} = csvread(strcat(outdomain, '-scalComp-model,', num2str(model_idx), '-rmse_std.csv'));
%         end
%     end
    

