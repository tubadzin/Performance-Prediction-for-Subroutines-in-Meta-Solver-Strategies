function model = learnModel(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names, varargin)
    if nargin < 7
        names = {};
    end

    if isfield(options, 'opt') && options.opt && size(X,1) > 1
        %== Set up inputs for DIRECT.
        opts.maxevals = options.hyp_opt_steps;
        opts.maxits = 1e10; % deactivated
        opts.showits = 1;
        bounds = repmat([options.paramsLowerBound, options.paramsUpperBound], length(options.tuning_params),1);

        %== Set up objective function for DIRECT.
        k = 2;
        func_handle = @(outer_model_params) evalOuterParamsCV(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names, k, outer_model_params, varargin{:});
        Problem.f = func_handle;

        %== Run DIRECT.
        [obj, optimized_outer_params, hist] = Direct(Problem, bounds, opts);
        optimized_outer_params
        obj
        options.tuning_params = optimized_outer_params;
    end
    model = learnModelWithFixedOuterParams(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names);
end

function model = learnModelWithFixedOuterParams(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names, varargin)
    %=== Wrapper around initializing a model and then preparing it for prediction.
    model = initModel(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names, varargin{:});
    if length(X) == 0
        if length(y) == 0
            model.y = 0.005; % minimal prediction
        end
        model.constant = true;
        model.options = options;
        return;
    end
    model = internalOptimizeModel(model, varargin{:});
    model = prepareModel(model);
end

function cv_performance = evalOuterParamsCV(X, y, cens, cat, catDomains, isclean, algo_deterministic, options, names, k, outer_model_params, varargin)
    %=== Combine performance of several folds calling evalOuterParams
    cv_performance = 0;
    options.tuning_params = outer_model_params;

    %== Use the same seed here all the time, but don't affect outer seed.
    s = rand('twister');
    rand('twister',1234);
    N = length(y);
    indices_for_cv = 1:N;
    randInd = randperm(N);
    rand('twister',s);

    %== Do cross-validation.
    startIdx = 1;
    for i=1:k
        endIdx = ceil(i*N/k);
        testIdx = startIdx:endIdx;
        trainIdx = setdiff(1:N, testIdx);
    
        cv_performance = cv_performance + evalOuterParams(X(trainIdx,:), y(trainIdx,:), cens(trainIdx,:), X(testIdx,:), y(testIdx,:), cens(testIdx,:), cat, catDomains, isclean, algo_deterministic, options, names, varargin{:});
        startIdx = endIdx+1;
    end
    cv_performance  = cv_performance/k;
end    

function performance = evalOuterParams(Xtrain, ytrain, censtrain, Xtest, ytest, censtest, cat, catDomains, isclean, algo_deterministic, options, names, varargin)
    %== Learn the model with the given outer params.
    model = learnModelWithFixedOuterParams(Xtrain, ytrain, censtrain, cat, catDomains, isclean, algo_deterministic, options, names, varargin{:});

    %== Evaluate the learned model.
    [yPredMean, yPredVar] = applyModel(model, Xtest, isclean);
    if options.logModel
        ytest = log10(ytest);
    end
    [rmse, ll, cc] = measures_of_fit(ytest, yPredMean, yPredVar, censtest);
    performance = rmse;
end