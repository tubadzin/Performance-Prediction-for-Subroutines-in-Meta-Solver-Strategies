function do_single_pass_of_var_imp(options_vec, savefilename, X, y, cat, catDomains, featureNames)
    perm = randperm(size(X,1));
    trainsize = ceil(length(perm)/5*3);
    validsize = ceil(length(perm)/5*1);
    train_idxs = perm(1:trainsize);
    valid_idxs = perm(trainsize+1:trainsize+validsize);
    test_idxs = perm(trainsize+validsize+1:end);
    fulltrain_idxs = [train_idxs,valid_idxs];
    for i=1:length(options_vec)
        options = options_vec{i};
        vimp.full_rmse(i) = get_rmse(X(fulltrain_idxs,:), y(fulltrain_idxs), X(test_idxs,:), y(test_idxs), cat, catDomains, options, featureNames, 1:size(X,2));
        
        numFeaturesToSelect = 10; %size(X,2);
        [selected_feature_idxs{i}, selected_feature_indices_in_order_added{i}, orig_rmses{i}] = fwdSelection(X(train_idxs,:), y(train_idxs), X(valid_idxs,:), y(valid_idxs), options, cat, catDomains, featureNames, numFeaturesToSelect);
        rmses{i} = [];
        for j=1:length(selected_feature_indices_in_order_added{i})
            sel_feats = sort(selected_feature_indices_in_order_added{i}(1:j));
            rmses{i}(end+1) = get_rmse(X(fulltrain_idxs,:), y(fulltrain_idxs), X(test_idxs,:), y(test_idxs), cat, catDomains, options, featureNames, sel_feats);
        end
        rmses{i}
        [tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatureNames] = use_feat_subset(X(train_idxs,:),X(test_idxs,:),cat,catDomains,featureNames,selected_feature_idxs{i});
        importance{i} = varImportance(tmpXTrain, y(train_idxs), tmpXValid, y(test_idxs), options, tmpCat, tmpCatDomains, tmpFeatureNames); %#ok<AGROW>
    end
    
    vimp.featureNames = featureNames;
    vimp.selected_feature_idxs = selected_feature_idxs;
    vimp.selected_feature_indices_in_order_added = selected_feature_indices_in_order_added;
    vimp.rmses = rmses;
    vimp.orig_rmses = orig_rmses;
    vimp.importance = importance;
    vimp.options_vec = options_vec;
	save(savefilename, 'vimp');
end

% function do_var_imp(options_vec, savefilename, X, y, cat, catDomains, featureNames)
%     perm = randperm(size(X,1));
%     trainsize = ceil(length(perm)/3*2);
%     train_idxs = perm(1:trainsize);
%     test_idxs = perm(trainsize+1:end);
%     for i=1:length(options_vec)
%         options = options_vec{i};
%         vimp.full_rmse(i) = get_rmse(X(train_idxs,:), y(train_idxs), X(test_idxs,:), y(test_idxs), cat, catDomains, options, featureNames, 1:size(X,2));
%         
%         numFeaturesToSelect = 10;
%         [selected_feature_idxs{i}, selected_feature_indices_in_order_added{i}, rmses{i}] = fwdSelection(X(train_idxs,:), y(train_idxs), X(test_idxs,:), y(test_idxs), options, cat, catDomains, featureNames, numFeaturesToSelect);
%         [tmpXTrain,tmpXValid,tmpCat,tmpCatDomains,tmpFeatureNames] = use_feat_subset(X(train_idxs,:),X(test_idxs,:),cat,catDomains,featureNames,selected_feature_idxs{i});
%         importance{i} = varImportance(tmpXTrain, y(train_idxs), tmpXValid, y(test_idxs), options, tmpCat, tmpCatDomains, tmpFeatureNames) %#ok<AGROW>
%     end
%     
%     vimp.featureNames = featureNames;
%     vimp.selected_feature_idxs = selected_feature_idxs;
%     vimp.selected_feature_indices_in_order_added = selected_feature_indices_in_order_added;
%     vimp.rmses = rmses;
%     vimp.importance = importance;
%     vimp.options_vec = options_vec;
% 	save(savefilename, 'vimp');
% end