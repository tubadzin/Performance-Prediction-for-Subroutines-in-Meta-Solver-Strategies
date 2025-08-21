function [yPred, yPredVar, samples_of_min] = applyMarginalModel(model, thetaTest, instanceFeatures, isclean, observationNoise, joint)

% if nargin < 3
%     filename = 'tmp_debug_applyMarginalModel_input.mat';
%     save(filename, 'model', 'thetaTest');
% elseif nargin < 4
%     filename = 'tmp_debug_applyMarginalModel_input_inst.mat';
%     save(filename, 'model', 'thetaTest', 'instanceFeatures');
% elseif nargin < 5
%     filename = 'tmp_debug_applyMarginalModel_input_inst_isclean.mat';
%     save(filename, 'model', 'thetaTest', 'instanceFeatures', 'isclean');
% elseif nargin < 6
%     filename = 'tmp_debug_applyMarginalModel_input_inst_isclean_obsNoise.mat';
%     save(filename, 'model', 'thetaTest', 'instanceFeatures', 'isclean', 'observationNoise');
% else
%     filename = 'tmp_debug_applyMarginalModel_input.mat';
%     save(filename, 'model', 'thetaTest', 'instanceFeatures', 'isclean', 'observationNoise', 'joint');
% end 

if isfield(model,'constant')
    if model.options.logModel
        yPred = ones(size(X,1),1) * log10(mean(model.y));
    else        
        yPred = ones(size(X,1),1) * mean(model.y);
    end
    yPredVar = ones(size(thetaTest,1),1) * var(model.y);
    samples_of_min = ones(size(thetaTest,1),1) * model.y(ceil(rand*length(model.y)));
    return
end
if nargin < 4
    isclean = 0;
end
if nargin < 5
    observationNoise = 0;
end
if nargin < 6
    joint = 0;
end
if ~isclean
    thetaTest = thetaTest(:,model.kept_Theta_columns);
end
if nargin < 3
    instanceFeatures = model.pcaed_instance_features;
    isclean = 1;
end
if ~isclean
    instanceFeatures = instanceFeatures(:,model.kept_X_columns);
    instanceFeatures = pca_fwd(instanceFeatures, model.sub, model.means, model.stds, model.pcVec, model.num_pca);
    isclean = 1;
end

if model.moreThanOneInstance == 0
    if nargout <= 2
        [yPred, yPredVar] = applyModel(model, thetaTest, isclean, observationNoise, joint);
    else
        [yPred, yPredVar, samples_of_min] = applyModel(model, thetaTest, isclean, observationNoise, joint);
    end
    return;
end
    
%=== Multiple instances => only random forest implemented.
if ~strcmp(model.type, 'rf')
    error('Multiple instances => only random forest implemented.');
end

[yPred, yPredVar, treemeans] = compute_from_leaves(model.module, thetaTest, instanceFeatures, model.algo_deterministic);

%=== Simple check to validate correctness by comparing to average of single
%predictions.
% for i=1:size(thetaTest,1)
%     [yPredSingle] = applyModel(model, [repmat(thetaTest(i,:),[2000,1]), instanceFeatures], isclean, observationNoise, joint);
%     yPredValid = mean(yPredSingle);
%     assert(abs(yPredValid - yPred(i)) < 1e-6);
% end


if model.options.logModel
    logtreemeans = log10(treemeans);
    yPred = mean(logtreemeans,2);
    yPredVar = var(logtreemeans,0,2);
    samples_of_min = min(logtreemeans,[],1);
else
    yPred = mean(treemeans,2);
    yPredVar = var(treemeans,0,2);
    samples_of_min = min(treemeans,[],1);
end
