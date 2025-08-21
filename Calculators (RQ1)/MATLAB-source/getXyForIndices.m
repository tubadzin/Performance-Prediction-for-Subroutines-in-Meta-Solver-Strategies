function [X,y] = getXyForIndices(Theta, all_features, matrix, theta_idxs, inst_idxs)
%=== Build X and y with light memory footprint.
X = [Theta(theta_idxs,:), all_features(inst_idxs,:)];
y = -ones(length(theta_idxs),1);
for i=1:length(theta_idxs)
    y(i) = matrix(theta_idxs(i), inst_idxs(i));
end
