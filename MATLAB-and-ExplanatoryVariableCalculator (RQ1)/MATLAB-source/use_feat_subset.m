function [tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatNames] = use_feat_subset(xTrain,xValid,cat,catDomains,featNames,selected_feature_idxs)
tmpCat = [];
tmpCatDomains = {};
for i=1:length(selected_feature_idxs)
    if find(cat==selected_feature_idxs(i))
        tmpCat = [tmpCat, i];
        tmpCatDomains{i} = catDomains{selected_feature_idxs(i)};
    end
end
tmpXTrain = xTrain(:,selected_feature_idxs);
tmpXValid = xValid(:,selected_feature_idxs);
tmpFeatNames = featNames(selected_feature_idxs);