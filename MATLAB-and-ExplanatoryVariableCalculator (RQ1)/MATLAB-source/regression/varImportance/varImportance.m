function importance=varImportance(xTrain, yTrain, xValid, yValid, options, cat, catDomains, featureNames)

%Do variable importance, by collecting stats for each leave-one out variable model.
n=size(xTrain,2);

fullRMSE = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featureNames, 1:n);

for i=1:n
    tmpRMSE = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featureNames, setdiff(1:n, i));
    imp(i) = tmpRMSE-fullRMSE; %#ok<NASGU>
    fprintf(strcat(['Without feature ', featureNames{i}, ': ', num2str(tmpRMSE), ' vs. ', num2str(fullRMSE), '(full)\n']));
end
imp = max(imp,0);
%importance = 100*imp/fullRMSE;
% importance = 100 * (imp ./ sum(imp));
if max(imp) == 0
    importance = imp;
else
    importance = 100*imp/max(imp);
end

assert(~any(isnan(importance)));
assert(~any(isinf(importance)));

% function rmse = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featNames, selected_feature_idxs)
% cens = zeros(length(yTrain),1);
% [tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatNames] = use_feat_subset(xTrain,xValid,cat,catDomains,featNames,selected_feature_idxs);    
% model = learnModel(tmpXTrain, yTrain, cens, tmpCat, tmpCatDomains, 0, 0, options, tmpFeatNames);
% [yPredValid] = applyModel(model, tmpXValid, 0, 0, 0);
% rmse = simple_measures_of_fit(options.logModel, yValid, yPredValid);
% 
% function [tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatNames] = use_feat_subset(xTrain,xValid,cat,catDomains,featNames,selected_feature_idxs)
% tmpCat = [];
% tmpCatDomains = {};
% for i=1:length(selected_feature_idxs)
%     if find(cat==selected_feature_idxs(i))
%         tmpCat = [tmpCat, i];
%         tmpCatDomains{i} = catDomains{selected_feature_idxs(i)};
%     end
% end
% tmpXTrain = xTrain(:,selected_feature_idxs);
% tmpXValid = xValid(:,selected_feature_idxs);
% tmpFeatNames = featNames(selected_feature_idxs);