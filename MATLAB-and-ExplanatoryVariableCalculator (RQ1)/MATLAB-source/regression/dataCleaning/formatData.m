function [xData, yData] = formatData(xData, yData, transformation)

if isfield(transformation, 'one_in_k_encoding') && transformation.one_in_k_encoding
	%=== Encode categorical parameters using a 1-in-k encoding, and then treat them as continuous afterwards.
	xData = one_in_k_encoding(xData, transformation.origCat, transformation.origCatDomains);
	cat = [];
	catDomains = {};
end

%=== Move the categorical columns first.
cont = setdiff(1:size(xData,2), transformation.cat);
tmpCat = xData(:,transformation.cat);
tmpCont = xData(:,cont);
xData = [tmpCat, tmpCont];

%=== Transform response variable.
yData = transformResponse(transformation.responseTransformation, yData);

%=== Special condition for models that do their own preprocessing/feature selection.
if isfield(transformation, 'only_one_in_k_and_y_trans') && transformation.only_one_in_k_and_y_trans
    return;
end

%=== Remove constant features.
xData = xData(:,transformation.kept_columns);

%=== Cat.parameters 
paramsCatxData=xData(:,[1:transformation.numCatParamsAfterRemovingConstants]);

%=== Split into parameters and instance features.
xData = xData(:,transformation.numCatParamsAfterRemovingConstants+1:end);


%=== PCA. (no PCA at this point)
% xData = xData*transformation.pcVec;


feat_select_by_dropping_linear_comb = (isfield(transformation, 'feat_predict_quadratic_kept'));
if feat_select_by_dropping_linear_comb 
    % Special case to conform to JACM paper protocol.
    keptfeatures = transformation.feat_predict_quadratic_kept;
    xDataDefault=xData*0+1;
    xDataDefault(find(xData==-512 | xData==-1024))=2;

    xDatafoo=[];
    xDefaultfoo = [];
    for i=1:length(keptfeatures)
        a =floor(keptfeatures(i)/100000);
        b =rem(keptfeatures(i), 100000);
        if a==0
            xDatafoo=[xDatafoo, xData(:,b)];
            xDefaultfoo = [xDefaultfoo, xDataDefault(:,b)];
        end
        if b==0
            xDatafoo=[xDatafoo, xData(:,a)];
            xDefaultfoo = [xDefaultfoo, xDataDefault(:,a)];
        end
        if a*b>0
            xDatafoo=[xDatafoo, xData(:,a).*xData(:,b)];
            xDefaultfoo = [xDefaultfoo, xDataDefault(:,a).*xDataDefault(:,b)]; % anything times broken is broken
        end
    end
    xData = xDatafoo;
    xDataDefault = xDefaultfoo;
    %=== Normalize data.
    xData = xData(:, transformation.nonconstant);
    if ~isempty(transformation.nonconstant)
        xData = (xData - repmat(transformation.bias, [size(xData,1),1])) ./ repmat(transformation.scale, [size(xData,1),1]);
        xDataDefault = xDataDefault(:, transformation.nonconstant);
        xData(find(xDataDefault>1))=0;
    end


else
    %=== Feature selection for linear instance features.
    foo=[paramsCatxData, xData];
    xData = foo(:,transformation.linearFeatureIndices);
    paramsCatxData = foo(:, transformation.linearCatParamsIndices);
    %=== deal with default feature
    xDataDefault=xData*0+1;
    xDataDefault(find(xData==-512 | xData==-1024))=2;
    %=== Build basis functions.
    %=== For data without parameters, this may include building quadratic features etc.
    %=== For data with parameters, it also builds combinations of parameters and instance features.
    [xData, NamesX, Aid, Bid] = buildBasisFunctions(xData, transformation.linearFeatureNames, transformation.doQuadratic);
    [xDataDefault, NamesX, Aid, Bid] = buildBasisFunctions(xDataDefault, transformation.linearFeatureNames, transformation.doQuadratic);

    %=== Normalize data.
    if ~isempty(transformation.nonconstant)
        xData = xData(:, transformation.nonconstant);
        xData = (xData - repmat(transformation.bias, [size(xData,1),1])) ./ repmat(transformation.scale, [size(xData,1),1]);
        xDataDefault = xDataDefault(:, transformation.nonconstant);
        xData(find(xDataDefault>1))=0;
    end

    finalxData = [paramsCatxData, xData];
    %=== Final feature selection to build small model.
    % finalFeatureIndices = transformation.setOfFinalFeatureIndices{transformation.bestIndexForFinalFeatures};
    % xData = finalxData(:,finalFeatureIndices);

    %=== Put categorical parameters first.
    tmpCatxTrain = finalxData(:,transformation.finalCatParamsIndices);
    tmpxTrain = finalxData(:,transformation.finalFeatureIndices);
    xData=[tmpCatxTrain, tmpxTrain];

    t = any(any(isnan(xData),2));
    t = t | any(any(isinf(xData),2));
    t = t | any(isnan(yData));
    t = t | any(isinf(yData));
    if any(t)
        error 'Empty values and infinity not allowed in either X, y.'
    end
end

if size(xData,2) == 0
    xData = zeros(size(xData,1),1);
end
