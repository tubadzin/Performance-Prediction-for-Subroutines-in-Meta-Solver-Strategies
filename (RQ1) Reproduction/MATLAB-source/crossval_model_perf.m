function [rmses, lls, ccs, cc_ranks, learnTimes, predTimes, y, all_y_cross, all_y_cross_var, cens] = crossval_model_perf(X, y, k, options, cat, catDomains, algo_deterministic, doplot, domain, axmin, axmax)
% 8 inputs ok
if nargin < 10
    axmin = 0.001;
    axmax = 3600;
    if nargin < 9
        domain = '';
        if nargin < 8
            doplot = 0;
        end
    end
end
%=== Do k-fold cross-validation to evaluate performance of different models.
N = length(y);
rand('twister',1234);
perm = randperm(N);
y=y(perm);
X=X(perm,:);

%=== Set up for simple model.
cens = zeros(N,1);               % counting censored data at the censoring threshold.
if options.logModel
    y = max(y, 0.005);
end
featureNames = {};

startIdx = 1;
rmses = [];
lls = [];
ccs = [];
cc_ranks = [];
learnTimes = [];
predTimes = [];
for i=1:k
    fprintf(strcat(['   Cross-validation ', num2str(i), '/', num2str(k), '... ']));
    endIdx = ceil(i*N/k);
    testIdx = startIdx:endIdx;
    trainIdx = setdiff(1:N, startIdx:endIdx);

    tic;
    model = learnModel(X(trainIdx,:), y(trainIdx), cens(trainIdx), cat, catDomains, 0, algo_deterministic, options, featureNames);
    learnTimes(i) = toc;

    tic;
    %=== Predict for test data of this fold.
    [y_cross, y_cross_var] = applyModel(model, X(testIdx,:), 0, 0, 0);
    all_y_cross(testIdx,:) = y_cross;
    all_y_cross_var(testIdx,:) = y_cross_var;
    predTimes(i) = toc;
    
    ytest = y(testIdx);
    censtest = cens(testIdx);
    if model.options.logModel
        [rmses(i), lls(i), ccs(i), cc_ranks(i)] = measures_of_fit(log10(ytest), y_cross, y_cross_var, censtest);
    else
        [rmses(i), lls(i), ccs(i), cc_ranks(i)] = measures_of_fit(ytest, y_cross, y_cross_var, censtest);
    end
%    [model.params(1),model.params(2),rmses(i)]

    startIdx = endIdx + 1;
end

if doplot
    if model.options.logModel
        [rmse, ll, cc] = measures_of_fit(log10(y), all_y_cross, all_y_cross_var, cens);
    else
        [rmse, ll, cc] = measures_of_fit(y, all_y_cross, all_y_cross_var, cens);
    end

    figure_prefix = strcat(['results/ehm/', domain, '_', options.unique_model_name]);
    title_prefix = '';
    plot_simple_pred_scatter(y, all_y_cross, all_y_cross_var, cens, rmse, cc, ll, figure_prefix, title_prefix, model.options.logModel, axmin, axmax)
end