function [rmse, ll, cc, learnTime, predTime] = simple_model_perf_train_test(Xtrain, ytrain, Xtest, ytest, cat, catDomains, options)

%=== Set up for simple model.
N = length(ytest);
censtrain = zeros(length(ytrain),1);               % counting censored data at the censoring threshold.
censtest = zeros(length(ytest),1);
if options.logModel
    ytrain = max(ytrain, 0.005);
    ytest = max(ytest, 0.005);
end
featureNames = {};

%=== Set up training and test indices.
rand('twister',1234);

tic;
model = learnModel(Xtrain, ytrain, censtrain, cat, catDomains, 0, options, featureNames);
learnTime = toc;

tic;
[y_pred, y_predvar] = applyModel(model, Xtest, 0, 0, 0);
predTime = toc;

if options.logModel
    [rmse, ll, cc] = measures_of_fit(log10(ytest), y_pred, y_predvar, censtest);
else
    [rmse, ll, cc] = measures_of_fit(ytest, y_pred, y_predvar, censtest);
end