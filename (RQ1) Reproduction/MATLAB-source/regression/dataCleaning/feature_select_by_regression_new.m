function [keepfeatures, history] = feature_select_by_regression_new(features, expansion, tor)
% using all other featues to predict one feature, it the prediction is
% really good, drop it since it is useless and can be obtained by the
% linear combination of all other features
% only support LR for Kevin
numfeat=size(features,2);
numinst=size(features,1);
best = 1;
history.drop=[];
history.acc=[];
keepfeatures=[];
maxIteration=10;
dropindx=1;
if expansion==0
    options.modelType='LR';
    options.logModel=0;
    keepfeatures=[1:numfeat];
    changed = 1;
    while changed
        changed = 0;
        % for each feature, use all others to predict the one
        tmpaccuracy=[];
        working_features = keepfeatures;
        for i = 1:length(keepfeatures)
            yid=keepfeatures(i);
            xid=setdiff(working_features, yid);
            % now build LR models
            model = learnModel(features(:, xid), features(:, yid), features(:, yid)*0, [], [], 1, 0, options, []); % do not do clean data
            % now make prediction
            predict = applyModel(model, features(:, xid), 1, [], []);
            % compute prediction accuracy on training data
            serror=sum((features(:, yid)-predict).^2);
            stotal=sum((features(:, yid)-mean( features(:, yid))).^2);
            rsquare =  1-serror/stotal;
            adjrsquare = 1 - (1-rsquare)*((length(features)-1)/(length(features)-length(xid)-1));
            % now recorde accuracy
            %             tmpaccuracy=[tmpaccuracy, adjrsquare];
            
            if adjrsquare > tor
                working_features = setdiff(working_features, yid);
                changed = 1;
                fprintf('Removing feature %d ...\n',yid);
            end
        end
        keepfeatures = working_features;
        %         best = max(tmpaccuracy);
        %         if best > tor
        %             bestid=find(tmpaccuracy==best);
        %
        %             history.drop=[history.drop, keepfeatures(bestid(1))];
        %             history.acc=[history.acc, best];
        %             keepfeatures=setdiff(keepfeatures, keepfeatures(bestid(1)));
        %         end
    end
end
if expansion==1
    options.modelType='LR';
    options.logModel=0;
    keepfeatures=zeros((numfeat+numfeat+numfeat*(numfeat-1)/2+1),1); % list of all possible features after expansion
    idx=1;
    
    for ii=1:numfeat+1
        for jj=ii:numfeat+1
            keepfeatures(idx)=(ii-1)*100000+(jj-1);
            idx=idx+1;
        end
    end
    keepfeatures=keepfeatures(2:end);
    predkeepfeatures=keepfeatures; % current allowed features for prediction

    %% now lets do forward selection
    for i=1:length(keepfeatures)
        %         fprintf('Try to drop (expanded) feature %d ...\n', keepfeatures(i));
        % Debug code; remove later.
        a =floor(keepfeatures(i)/100000);
        b =rem(keepfeatures(i), 100000);

        if a==0
            onefeature=features(:,b);
        end
        if b==0
            onefeature=features(:,a);
        end
        if a*b>0
            onefeature=features(:,a).*features(:,b);
        end
        otherfeatureid=setdiff(predkeepfeatures,keepfeatures(i));

        X = ones(length(onefeature),1);
        delta=1e-2; % = 1e-16;
        A = X'*X+delta;
        Ainv = inv(A);
        XtY = X' * onefeature;
        done=0;
        nowfeatures = [];
        %% Forward selection up to maxIteration many predictors
        for t=1:maxIteration
            bestacc=inf*-1;
            result=[];
            if done>0
                continue;
            end
            for tt=1:length(otherfeatureid) %- size(features,2)
                if (~isempty(find(nowfeatures==otherfeatureid(tt))))
                    continue;
                end
                if done>0
                    continue;
                end
                nowfeatures(t)=otherfeatureid(tt);
                fooa=floor(nowfeatures(t)/100000);
                foob=rem(nowfeatures(t), 100000);

                if fooa==0
                    foofeat=features(:,foob);
                end
                if foob==0
                    foofeat=features(:,fooa);
                end
                if fooa*foob>0
                    foofeat=features(:,fooa).*features(:,foob);
                end

                % construct the model with this feature added
                model = modelAddFeature(X, foofeat, onefeature, A, Ainv, delta, XtY);

                % record the RMSE and MAE on validation data for this feature
                newX = [X foofeat];
                preds = newX*model;

                serror=sum((onefeature-preds).^2);
                stotal=sum((onefeature-mean( onefeature)).^2);
                rsquare =  1-serror/stotal;
                adjrsquare = 1 - (1-rsquare)*(numinst-1)/(numinst-t-1);
                result(tt)=adjrsquare;

                bestacc = max(result(tt));
                if (bestacc == adjrsquare)
                    best_subset = otherfeatureid(tt);
                end
                if bestacc >=tor
                    Afoo = newX' * newX + delta * eye(size(newX,2));
                    Ainvfoo=pinv(Afoo);
                    XtYfoo = newX' * onefeature;
                    modelfoo=Ainvfoo * XtYfoo;
                    preds = newX*modelfoo;
                    serror=sum((onefeature-preds).^2);
                    stotal=sum((onefeature-mean(onefeature)).^2);
                    rsquare =  1-serror/stotal;
                    adjrsquare = 1 - (1-rsquare)*(numinst-1)/(numinst-t-1);
                    if adjrsquare > tor
                        done=1;
                        fprintf('Removing feature %d with cc %f ...\n',[keepfeatures(i), adjrsquare]);
                        history.drop(dropindx)=keepfeatures(i);
                        history.acc(dropindx)=adjrsquare;
                        dropindx=dropindx+1;
                        break;
                    end
                end

            end
            fooa=floor(best_subset/100000);
            foob=rem(best_subset, 100000);
            if fooa==0
                foofeat=features(:,foob);
            end
            if foob==0
                foofeat=features(:,fooa);
            end
            if fooa*foob>0
                foofeat=features(:,fooa).*features(:,foob);
            end
            X = [X foofeat];

            %%% replace those line with 3 lines up
            A = X' * X + delta * eye(size(X,2));
            Ainv=pinv(A);
            XtY = X' * onefeature;
            nowfeatures(t)=best_subset;
            %             model=Ainv * XtY;
            %             preds = X*model;
            %             % compute prediction accuracy on training data
            %             serror=sum((onefeature-predict).^2);
            %             stotal=sum((onefeature-mean(onefeature)).^2);
            %             rsquare =  1-serror/stotal;
            %             adjrsquare = 1 - (1-rsquare)*(numinst-1)/(numinst-t));
            %             bestacc=adjrsquare;
            %
        end
        if done ==1
            predkeepfeatures=otherfeatureid;
        end
    end

    keepfeatures=predkeepfeatures;
end