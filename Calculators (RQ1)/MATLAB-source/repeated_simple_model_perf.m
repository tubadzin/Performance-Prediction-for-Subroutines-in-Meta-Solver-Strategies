function [rmses, lls, ccs, cc_ranks, learnTimes, predTimes] = repeated_simple_model_perf(X, y, numTrain, options, cat, catDomains, algo_deterministic, numRepeats)
numTrain = min(numTrain, ceil(length(y)*9/10.0));
rmses = -ones(numRepeats,1);
lls = -ones(numRepeats,1);
ccs = -ones(numRepeats,1);
learnTimes = -ones(numRepeats,1);
predTimes = -ones(numRepeats,1);

%=== Set up for simple model.
N = length(y);
cens = zeros(N,1);               % counting censored data at the censoring threshold.
if options.logModel
    y = max(y, 0.005);
end
featureNames = {};

for i=1:numRepeats
    rand('twister',i);
    
    %=== Set up training and test indices.
    perm=randperm(N);
    trainIdx = perm(1:numTrain);
    testIdx = perm(numTrain+1:end);
%    testIdx = testIdx(1:200);

    tic;
    model = learnModel(X(trainIdx,:), y(trainIdx), cens(trainIdx), cat, catDomains, 0, algo_deterministic, options, featureNames);
    learnTimes(i) = toc;

    tic;
    [y_pred, y_predvar] = applyModel(model, X(testIdx,:), 0, 0, 0);
    predTimes(i) = toc;

    ytest = y(testIdx);
    censtest = cens(testIdx);
    if options.logModel
        [rmses(i), lls(i), ccs(i), cc_ranks(i)] = measures_of_fit(log10(ytest), y_pred, y_predvar, censtest);
    else
        [rmses(i), lls(i), ccs(i), cc_ranks(i)] = measures_of_fit(ytest, y_pred, y_predvar, censtest);
    end
end