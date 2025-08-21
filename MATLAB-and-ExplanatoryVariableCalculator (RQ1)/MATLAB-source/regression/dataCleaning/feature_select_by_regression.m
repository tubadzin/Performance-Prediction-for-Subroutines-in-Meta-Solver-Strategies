function [kept_features, history] = feature_select_by_regression(features, type, tor, max_to_keep)
% using all other featues to predict one feature, it the prediction is
% really good, drop it since it is useless and can be obtained by the
% linear combination of all other features
% only support LR for Kevin
if nargin < 4
    max_to_keep = inf; % de-activate this criterion
end
numfeat=size(features,2);
best = inf;
history.drop=[];
history.acc=[];
kept_features=[1:numfeat];
if strcmp(type, 'LR')
    options.modelType=type;
    options.logModel=0;    
    count = 1;
    while best > tor || length(kept_features) > max_to_keep
        % for each feature, use all others to predict the one
        tmpaccuracy=[];
        for i = 1:length(kept_features)
            yid=kept_features(i);
            xid=setdiff(kept_features, yid);
            % now build LR models 
            model = learnModel(features(:, xid), features(:, yid), features(:, yid)*0, [], [], 1, 0, options, []); % do not do clean data
            % now make prediction 
            predict = applyModel(model, features(:, xid), 1, [], []);
            % compute prediction accuracy on training data
%             serror=sum((features(:, yid)-predict).^2);
%             stotal=sum((features(:, yid)-mean( features(:, yid))).^2);
%             rsquare =  1-serror/stotal;
%             adjrsquare = 1 - (1-rsquare)*((length(features)-1)/(length(features)-length(xid)));
            adjrsquare = corr(features(:, yid), predict);
            if isnan(adjrsquare)
                adjrsquare = 0;
            end
            % now recorde accuracy
            tmpaccuracy=[tmpaccuracy, adjrsquare];            
        end
        best = max(tmpaccuracy);
        if best > tor || length(kept_features) > max_to_keep
            bestid=find(tmpaccuracy==best);
            
            history.drop=[history.drop, kept_features(bestid(1))];
            history.acc=[history.acc, best];
            kept_features=setdiff(kept_features, kept_features(bestid(1)));
            fprintf(['Iteration ', num2str(count), ': dropping feature with rho = ', num2str(best), '\n']);
%                 ' with current index ', num2str(bestid(1)), '\n']);
        end
        count = count+1;
    end
else
    fprintf('Error, %s is not supported \n', type);
end
kept_features