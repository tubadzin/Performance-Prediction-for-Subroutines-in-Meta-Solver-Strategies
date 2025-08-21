function [yPred, yPredVar, timeToPredict] = modelPredictForJava(model, Theta, InstanceFeatures)

tic;
X = [Theta, InstanceFeatures];
[yPred, yPredVar] = applyModel(model, X, 0, 0, 0);
timeToPredict = toc;