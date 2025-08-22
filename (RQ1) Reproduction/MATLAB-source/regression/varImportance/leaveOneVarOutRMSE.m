function rmse=leaveOneVarOutRMSE(xTrain, yTrain, xValid, yValid, options, cat, catDomains, featureNames)

n = size(xTrain,2);
for i=1:n
    rmse(i,1) = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featureNames, setdiff(1:n, i));
end
rmse(n+1,1) = get_rmse(xTrain, yTrain, xValid, yValid, cat, catDomains, options, featureNames, 1:n);
rmse