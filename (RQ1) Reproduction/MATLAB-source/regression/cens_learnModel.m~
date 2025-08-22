% Theta_inst_idxs is a N*2 matrix of [theta_idx, inst_idx] pairs indexing into ThetaUniqSoFar and all_instance_features.
function model = cens_learnModel(Theta_inst_idxs, ThetaUniqSoFar, all_instance_features, y, cens, thetaCat, thetaCatDomains, xCat, xCatDomains, cond_params_idxs, parent_param_idxs, ok_parent_value_idxs, isclean, options, names, varargin)
if isfield(options, 'censoring') && (strcmp(options.censoring, 'fill_in_mean') || strcmp(options.censoring, 'fill_in_samples'))
    %=== Get indices of censored and noncensored data points.
    noncens_idx = find(cens==0);
    
    if isempty(noncens_idx) % if all censored predict max. value
        model.constant = true;
        model.options = options;
        if options.logModel
            model.y = ones(length(y),1) * log10(model.options.kappa_max * model.options.cutoff_penalty_factor);
        else
            model.y = ones(length(y),1) * model.options.kappa_max * model.options.cutoff_penalty_factor;
        end
        return;
    end
    
    %=== Learn model on noncensored data.
    model = cens_internalLearnModel(Theta_inst_idxs(noncens_idx,:), ThetaUniqSoFar, all_instance_features, y(noncens_idx), cens(noncens_idx), thetaCat, thetaCatDomains, xCat, xCatDomains, cond_params_idxs, parent_param_idxs, ok_parent_value_idxs, isclean, options, names, [], [], varargin);

    %=== Set up variables that don't change during the data imputation.
    cens_idx = find(cens==1);
    XtrainCens = [ThetaUniqSoFar(Theta_inst_idxs(cens_idx,1),:), all_instance_features(Theta_inst_idxs(cens_idx,2),:)];
    
    seed_to_learn_and_impute = ceil(rand*10000000000);
    %=== Iteratively get predictive distribution and re-fit model.
    old_prng_state = rand('twister');
    for imputation_iteration = 1:5
        rand('twister', seed_to_learn_and_impute); % same seed in each iteration to reduce variance.

        %== Get predictive distribution.
        [censMeanPred, censVarPred] = applyModel(model, XtrainCens, 0, 0, 0);
        if imputation_iteration == 1
            old_mean_imputed = censMeanPred;
        end
        
        %== Re-fit.
        model = cens_internalLearnModel(Theta_inst_idxs, ThetaUniqSoFar, all_instance_features, y, cens, thetaCat, thetaCatDomains, xCat, xCatDomains, cond_params_idxs, parent_param_idxs, ok_parent_value_idxs, isclean, options, names, censMeanPred, censVarPred, varargin);

        %== Check how much the imputed values changed; break if too little 
%             average_increase_in_imputed_values = mean(censMeanPred - oldCensMeanPred);
%             fprintf(['Average increase of imputed values in imputation_iteration ', num2str(imputation_iteration-1), ':', num2str(average_increase_in_imputed_values), '\n']);
%         average_increase_in_imputed_values = mean(model.mean_imputed - censMeanPred);
%         oldCensMeanPred = censMeanPred;
%         fprintf(['Alternative computation of average increase of imputed values in imputation_iteration ', num2str(imputation_iteration), ':', num2str(average_increase_in_imputed_values), '\n']);
        average_increase_in_imputed_values = mean(model.mean_imputed-old_mean_imputed);
        fprintf(['Mean of imputed values in imputation_iteration ', num2str(imputation_iteration), ':', num2str(average_increase_in_imputed_values), '\n']);
        old_mean_imputed = model.mean_imputed;
        if average_increase_in_imputed_values < 1e-10 && imputation_iteration >= 2
            fprintf(['Mean of imputed values stopped increasing in imputation_iteration ', num2str(imputation_iteration), '.']);
            return;
        end
    end
    rand('twister', old_prng_state);
else
    model = cens_internalLearnModel(Theta_inst_idxs, ThetaUniqSoFar, all_instance_features, y, cens, thetaCat, thetaCatDomains, xCat, xCatDomains, cond_params_idxs, parent_param_idxs, ok_parent_value_idxs, isclean, options, names, [], [], varargin);
end