function [yPred, yPredVar] = cens_applyModel(model, X, isclean, observationNoise, joint)
if nargin < 5
    joint = 0;
end
if nargin < 4
    observationNoise = 0;
end
if nargin < 3
    isclean = 0;
end

if isfield(model,'constant')
    if model.options.logModel
        yPred = ones(size(X,1),1) * log10(mean(model.y));
    else        
        yPred = ones(size(X,1),1) * mean(model.y);
    end
    yPredVar = ones(size(X,1),1) * var(model.y);
    yPredVar = max(yPredVar, model.options.min_variance);
    return
end
if ~isclean
    Theta = X(:, 1:model.numParameters);
    instanceFeatures = X(:, model.numParameters+1:end);
    
    Theta = Theta(:, model.kept_Theta_columns);
    
    if isempty(model.kept_X_columns)
        instanceFeatures = zeros(size(instanceFeatures,1),1); % dummy instance feature
    else
        instanceFeatures = instanceFeatures(:, model.kept_X_columns);
    end
    
    instanceFeatures = pca_fwd(instanceFeatures, model.sub, model.means, model.stds, model.pcVec, model.num_pca);
    
    X = [Theta, instanceFeatures];
    isclean = 1;
end


if strcmp(model.type, 'rf')
    samples_of_min = inf * ones(1, length(model.module));
    numchunks = ceil(size(X,1)*length(model.module) / 1000000);
    lenchunks = ceil(size(X,1)/numchunks);
    lower = 1;
    for i=1:numchunks
        upper = min([lower+lenchunks, size(X,1)]);
        subX = X(lower:upper,:);

        %=== Regular prediction for subX.
        treemeans = zeros(size(subX,1), length(model.module));
        treevars  = zeros(size(subX,1), length(model.module));
        for m=1:length(model.module)
            [treemeans(:,m), treevars(:,m)] = fh_simple_one_treeval(model.module{m}.T, subX, model.options.strategyForMissing);
            assert(~ any(isnan(treemeans(:,m))));
            if any( (isinf(treemeans(:,m))) )
                subX
                sub_treemeans = treemeans(:,m)
                T_m = model.module{m}.T
                assert(~ any(isinf(treemeans(:,m))));
            end
        end
        subyPred = mean(treemeans,2);

        subyPredVar = mean(treemeans.^2 + treevars, 2) - subyPred.^2;
        if any(subyPredVar < -1e-4)
            subyPredVar
            error('variance < -1e-4');
        end
        subyPredVar = max(subyPredVar, model.options.min_variance);
%             subyPredVar = max(subyPredVar, 1e-14);
%             subyPredVar = var(treemeans,0,2);
        subsamples_of_min = min(treemeans,[],1);

        %=== Regular prediction for subX.
        yPred(lower:upper,1) = subyPred;
        yPredVar(lower:upper,1) = subyPredVar;
        samples_of_min = min([samples_of_min;subsamples_of_min],[],1);
        lower = upper+1;
        if lower > size(X,1)
            break;
        end
    end

    assert(~any(isnan(yPredVar)));
    assert(~any(isinf(yPredVar)));
    return
end

if isempty(strfind(model.type, 'GP'))
    error ('No such model type defined yet!');
end

global gprParams;
gprParams = [];
gprParams.combinedcat = model.cat;
gprParams.combinedcont = model.cont;
%             gprParams.algoParam = model.algoParam;

%=== Normalize X using same normalization as before:
%             X = X(:, model.good_feats);
%             X = X - repmat(model.means, [size(X,1),1]);
%             X = X ./ repmat(model.stds, [size(X,1),1]);

if isfield(model.options, 'ppSize') && model.options.ppSize > 0
    if joint
        error 'Joint predictions are not implemented yet for SRPP.'
    end

    %=== GP with subset of regressors or projected process
    [yPred, S2SR, S2PP] = gprSRPPfwd(model.Kmm, model.invKmm, model.saved_1, model.saved_2, model.params, model.covfunc, model.pp_index, model.X, X, observationNoise);
    yPredVar = S2PP;
    if (any(yPredVar < 0))
        debug_filename = 'debug_file_for_neg_var_in_pp.mat';
        bout(sprintf(strcat(['\n\nWARNING: predicted variance is negative: ', num2str(min(yPredVar)), ', saving workspace to ', debug_filename])));
        save(debug_filename);
    end
    yPredVar = max(yPredVar, 1e-10);
%                yPredVar = S2SR; % the two seem very similar, but the GP book suggests PP is better far away from the data.
else
    %=== Normal GP
    if isfield(model, 'useCensoring')
        [yPred, yPredVar] = gprCensorFwd(model.X, model.L_nonoise, model.alpha_nonoise, model.invK_times_invH_times_invK, model.covfunc, model.params, X, observationNoise, joint);
    else
        [yPred, yPredVar] = gprFwd(model.X, model.L, model.invL, model.alpha, model.covfunc, model.params, X, observationNoise, joint);
    end
end

yPredVar = max(yPredVar, model.options.min_variance);
if isfield(model, 'meanToAdd')
    yPred = yPred + model.meanToAdd;
end