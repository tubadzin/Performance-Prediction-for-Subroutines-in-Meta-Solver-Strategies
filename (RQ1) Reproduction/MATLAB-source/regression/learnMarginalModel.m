function model = learnMarginalModel(Theta, used_instance_idxs, all_instance_features, y, cens, thetaCat, thetaDomains, xCat, xCatDomains, isclean, algo_deterministic, al_opts, names, varargin)
% If al_opts.logModel = 1, then we logify the data before fitting the model: y = log10(max(y,0.005))
if nargin < 9
    names = {};
end

% filename = 'tmp_debug_learnMarginalModel_input.mat';
% save(filename, 'Theta', 'used_instance_idxs', 'all_instance_features', 'y', 'cens', 'thetaCat', 'thetaDomains', 'xCat', 'xCatDomains', 'isclean', 'algo_deterministic', 'al_opts', 'names');

% %=== Diagnostic output
% Theta
% sizeOfTheta = size(Theta)
% 
% used_instance_idxs
% sizeOfused_instance_idxs = size(used_instance_idxs)
% 
% all_instance_features
% sizeOfall_instance_features = size(all_instance_features)
% 
% y
% sizeOfy = size(y)
% 
% cens
% thetaCat
% sizeOfthetaCat = size(thetaCat)
% 
% thetaDomains
% celldisp(thetaDomains)
% sizeOfthetaDomains = size(thetaDomains)
% 
% xCat
% sizeOfxCat = size(xCat)
% 
% xCatDomains
% celldisp(xCatDomains)
% sizeOfxCatDomains = size(xCatDomains)
% 
% isclean
% algo_deterministic
%    ------------------------------------------------------------//

if(~isempty(xCat))
    error('so far, this code only allows numerical features');
end

if ~isclean
    constant_Theta_columns = determine_constants(Theta, 1);
    kept_Theta_columns = setdiff(1:size(Theta,2), constant_Theta_columns);

    %=== Update thetaCat after dropping some columns of theta.
    %=== E.g. thetaCat = [3,7,8] and params 2&3 are constant. New thetaCat is [5,6]
    is_cat_theta = zeros(size(Theta,2),1);
    is_cat_theta(thetaCat) = 1;
    is_cat_theta = is_cat_theta(kept_Theta_columns);
    thetaCat = find(is_cat_theta);
    thetaDomains = thetaDomains(kept_Theta_columns);
    
    Theta = Theta(:,kept_Theta_columns);

    constant_X_columns = determine_constants(all_instance_features(used_instance_idxs,:), 1);
    kept_X_columns = setdiff(1:size(all_instance_features,2), constant_X_columns);
    all_instance_features = all_instance_features(:,kept_X_columns);
    
    numPCA = al_opts.pca;
    [pcaed_instance_features, sub, means, stds, pcVec] = do_pca(all_instance_features, numPCA);
    isclean = 1;
end

PiFeat = pcaed_instance_features(used_instance_idxs,:);

al_opts.marginalModel = 1; % only difference: numAlgoParam
%bout(sprintf(strcat(['Learning model with ', num2str(size(Theta,1)), ' data points of dimension ', num2str(size(Theta, 2) + size(PiFeat,2)),' ... \n'])));
%num_features = size(pcaed_instance_features,2) % mc: check that PCA performed on instance features, size is 7

model = learnModel([Theta,PiFeat], y, cens, [thetaCat, length(thetaCat) + xCat], [thetaDomains,  xCatDomains], isclean, algo_deterministic, al_opts, names, varargin(:));

model.moreThanOneInstance = (size(all_instance_features,1) > 1) && (size(all_instance_features,2) > 0);

model.kept_Theta_columns = kept_Theta_columns;
model.kept_X_columns = kept_X_columns;
model.sub = sub;
model.means = means;
model.stds = stds;
model.pcVec = pcVec;
model.num_pca = numPCA;

% model.trainInstanceFeatures = pca_fwd(all_instance_features, model.sub, model.means, model.stds, model.pcVec, model.num_pca);
model.pcaed_instance_features = pcaed_instance_features;
