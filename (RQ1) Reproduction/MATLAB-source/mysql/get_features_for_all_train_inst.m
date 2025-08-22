function all_features_for_training_instances = get_features_for_all_train_inst(func, ehFeature, tuningScenario)
% Get the features for all training instances. Features can be computed
% once in an offline stage, and then simply be read in.
% In contrast to per-instance approaches, in (static) algorithm
% configuration, we never need to compute features online.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ FEATURES FOR ALL TRAINING INSTANCES ONCE AND FOR ALL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
bout(sprintf(['Reading features for all instances from DB, once and for all...\n']));
N = func.numTrainingInstances;

if ehFeature < 4
    if func.cheap || func.matlab_fun
        all_features_for_training_instances = (1:N)';
        xNames = {'feat'};
    else
        assert(~func.matlab_fun);
        all_instIDs = getTrainID((1:N)', func);

        [all_features_for_training_instances, xNames] = getFeaturesForID(func.mip, all_instIDs);
    end
else
    all_features_for_training_instances = zeros(N,0);
    xNames = {};
end

if strcmp(tuningScenario, 'SPEAR-ibm-al-non-satelitesolved') || strcmp(tuningScenario, 'SPEAR-swv-al-non-satelitesolved')
    all_features_for_training_instances = all_features_for_training_instances(:, [1:44,46:end]);
end

if size(all_features_for_training_instances,1) == 1 && size(all_features_for_training_instances,2) > 0
    error 'If there is only one instance, instance features do not make any sense.';
end

assert(~any(any(isnan(all_features_for_training_instances))));
assert(~any(any(all_features_for_training_instances==-512)));

%=== Compute default feature.
yDef = [];
if ehFeature == 1 || ehFeature == 3
    s=rand('twister');
    %=== Get the default performance on each instance and use it as a feature.
    bout(sprintf(strcat(['\n\nGetting performance of the default on ', num2str(N), ' training instances.\n\n'])));
    theta_idx = update_if_new_param_config(func,func.default_values');
    defUsedSeeds = getTrainSeedsForInstanceNumbers(func, theta_idx*ones(N,1), (1:N)');
    [yDef, defCens] = func.funcHandle(theta_idx*ones(N,1), (1:N)', defUsedSeeds, func.cutoff * ones(N,1), tuningScenario);
    rand('twister', s);
    if ~(func.cheap || func.matlab_fun)
        yDef = max(yDef,0.005);
    end
end
if isempty(yDef)
    yDef = zeros(N,1);
end

switch ehFeature
    case 1
        all_features_for_training_instances = [all_features_for_training_instances, yDef];
    case 3
        all_features_for_training_instances = yDef;
    case 4
        all_features_for_training_instances = (1:N)';
end

readin_time = toc;
bout(sprintf(['Reading features for all instances from DB, once and for all...  took ', num2str(readin_time), 's.\n']));