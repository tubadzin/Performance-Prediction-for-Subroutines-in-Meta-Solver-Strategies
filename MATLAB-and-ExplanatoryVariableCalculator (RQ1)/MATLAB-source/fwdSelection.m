function [selected_feature_idxs,selected_feature_idxs_in_order_added,rmses]=fwdSelection(xTrain, yTrain, xValid, yValid, options, cat, catDomains, featureNames, numFeaturesToSelect)
rmses = [];

%Do forward feature selection in a model-independent way.
n=size(xTrain,2);
% if numFeaturesToSelect >= n
%     selected_feature_idxs = (1:n);
%     selected_feature_idxs_in_order_added = (1:n);
%     return;
% end

selected_feature_idxs = [];
%emptyModelRMSE = get_rmse(xTrain, yTrain, cens, cat, catDomains, options, featureNames, selected_feature_idxs);

selected_feature_idxs_in_order_added = [];
while length(selected_feature_idxs) < numFeaturesToSelect
    length(selected_feature_idxs)
    unpicked_features = setdiff(1:n, selected_feature_idxs);
    rmse = inf*ones(length(unpicked_features),1);
    for i=1:length(unpicked_features)
        var_idx = unpicked_features(i);
        
        %=== Don't waste time with constant features, set their rmse to inf
        %=== so they are not picked by fwd_sel.
%         x = [xTrain(:,var_idx); xValid(:,var_idx)];
        x = xTrain(:,var_idx); 
        if std(x) < 1e-5
            rmse(i) = inf;
        else
            rmse(i) = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featureNames, sort([selected_feature_idxs, var_idx]));
        end
    end
	min_rmse = min(rmse)
    min_idxs = find(rmse <= min_rmse + 1e-5);
    best_idx = min_idxs(ceil(rand*length(min_idxs)));
    selected_feature_idxs = [selected_feature_idxs, unpicked_features(best_idx)]; %#ok<AGROW>
    rmses(end+1) = min_rmse;
    selected_feature_idxs_in_order_added(end+1) = unpicked_features(best_idx);
    selected_feature_idxs = sort(selected_feature_idxs); % need to sort to make sure cat. features come before noncat., even if added later
end

