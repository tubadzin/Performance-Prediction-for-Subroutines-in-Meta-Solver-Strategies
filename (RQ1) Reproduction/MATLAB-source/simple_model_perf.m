function [rmse, ll, cc, learnTime, predTime] = simple_model_perf(X, y, numTrain, options)

%=== Set up for simple model.
N = length(y);
cens = zeros(N,1);               % counting censored data at the censoring threshold.
cat = [];                        % no categorical inputs
catDomains = [];                 % still none
if options.logModel
    y = max(y, 0.005);
end
featureNames = {};

%=== Set up training and test indices.
rand('twister',1234);
perm=randperm(N);
trainIdx = perm(1:numTrain);
testIdx = perm(numTrain+1:end);

tic;
model = learnModel(X(trainIdx,:), y(trainIdx), cens(trainIdx), cat, catDomains, 0, options, featureNames);
learnTime = toc;

tic;
[y_pred, y_predvar] = applyModel(model, X(testIdx,:), 0, 0, 0);
predTime = toc;

ytest = y(testIdx);
censtest = cens(testIdx);
if options.logModel
    [rmse, ll, cc] = measures_of_fit(log10(ytest), y_pred, y_predvar, censtest);
else
    [rmse, ll, cc] = measures_of_fit(ytest, y_pred, y_predvar, censtest);
end