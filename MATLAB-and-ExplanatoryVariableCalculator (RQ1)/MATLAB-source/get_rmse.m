function rmse = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featNames, selected_feature_idxs)
rand('twister',1234); % fix the random seed for variance reduction: e.g., in RFs, all feature subsets get the same random subsamples 
cens = zeros(length(yTrain),1);
% tmpXTrain = xTrain(:,selected_feature_idxs);
% tmpXValid = xValid(:,selected_feature_idxs);
% tmpFeatureNames = featureNames(selected_feature_idxs);
% 
% tmpCat = intersect(cat, selected_feature_idxs);
% idxs = [];
% for i=1:length(tmpCat)
%     idx = find(cat==tmpCat(i));
%     idxs = [idxs, idx]; %#ok<AGROW>
% end
% tmpCatDomains = catDomains(idxs);
% assert(length(tmpCat) == length(cat(idxs)));
% for i=1:length(tmpCat)
%     assert(tmpCat(i), cat(idxs(i)));
% end
[tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatNames] = use_feat_subset(xTrain,xValid,cat,catDomains,featNames,selected_feature_idxs);    

model = learnModel(tmpXTrain, yTrain, cens, tmpCat, tmpCatDomains, 0, 0, options, tmpFeatNames);
[yPredValid] = applyModel(model, tmpXValid, 0, 0, 0);
rmse = simple_measures_of_fit(options.logModel, yValid, yPredValid);