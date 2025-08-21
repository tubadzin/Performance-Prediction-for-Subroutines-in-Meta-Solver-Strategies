function vimp = cv_varImportance(X, y, options, cat, catDomains, featureNames)

rand('twister',1234);
N = length(y);
perm = randperm(N);
y=y(perm);
X=X(perm,:);

k=10;
startIdx = 1;
importance = zeros(size(X,2), k);
rmses = -ones(size(X,2)+1,k);
for i=1:k
    fprintf(strcat(['   Cross-validation for importance of input ', num2str(i), '/', num2str(k), '... ']));
    endIdx = ceil(i*N/k);
    testIdx = startIdx:endIdx;
    trainIdx = setdiff(1:N, startIdx:endIdx);

    numFeaturesToSelect = 10;
%     selected_feature_idxs = fwdSelection(X(trainIdx,:), y(trainIdx), X(testIdx,:), y(testIdx), options, cat, catDomains, featureNames, numFeaturesToSelect);
% %     cens = zeros(length(trainIdx),1);
% %     model = learnModel(X(trainIdx,selected_feature_idxs), y(trainIdx), cens, cat, catDomains, 0, 0, options, featureNames);
% %     [yPred] = applyModel(model, X(testIdx,selected_feature_idxs), 0, 0, 0);
% %     rmse = simple_measures_of_fit(options.logModel, y(testIdx), yPred)
%     importance(selected_feature_idxs,i) = varImportance(X(trainIdx,selected_feature_idxs), y(trainIdx), X(testIdx,selected_feature_idxs), y(testIdx), options, cat, catDomains, featureNames); %#ok<AGROW>

%     options.nSub = 100;
%     importance(1:size(X,2),i) = varImportance(X(trainIdx,:), y(trainIdx), X(testIdx,:), y(testIdx), options, cat, catDomains, featureNames) %#ok<AGROW>

%     options.nSub = 100;
    rmses(:,i) = leaveOneVarOutRMSE(X(trainIdx,:), y(trainIdx), X(testIdx,:), y(testIdx), options, cat, catDomains, featureNames)
    startIdx = endIdx + 1;
end
reference = repmat(rmses(end,:), [size(X,2),1]);
importance = (reference-rmses(1:end-1,:)) ./ reference;
importance = max(importance, 0);

mean_importance = mean(importance,2);
std_importance = std(importance,0,2);

[tmp, idx] = sort(mean_importance);
%k = min(20, length(idx));
vimp.idx=idx;
vimp.feature = featureNames(idx);
vimp.mean_importance = mean_importance(idx);
vimp.std_importance = std_importance(idx);
vimp.importance = importance(idx,:);
vimp.rmses = rmses;