function [rmse, ll, cc] = get_results_for_quadrant( model, Theta, all_features, matrix, nTest, cand_theta_idxs, cand_inst_idxs )
new_cand_theta_idxs = repmat(cand_theta_idxs, [length(cand_inst_idxs), 1]);
new_cand_theta_idxs = new_cand_theta_idxs(:);
new_cand_inst_idxs = repmat(cand_inst_idxs', [length(cand_theta_idxs), 1]);

perm = randperm(length(new_cand_theta_idxs));
actual_nTest = min(length(perm), nTest);
test_idx = perm(end-actual_nTest+1:end);

start_idx = 1;
tic;
chunk_size = 10000;
while start_idx <= length(test_idx)
    end_idx = min(start_idx+chunk_size, length(test_idx));
    this_idxs = test_idx(start_idx:end_idx);
    [Xtest,ytest(start_idx:end_idx,1)] = getXyForIndices(Theta, all_features, matrix, new_cand_theta_idxs(this_idxs), new_cand_inst_idxs(this_idxs));
    [ypred(start_idx:end_idx,1), ypredvar(start_idx:end_idx,1)] = applyModel(model, Xtest, 0, 0, 0);
    start_idx = end_idx+1
end

% [Xtest,ytest] = getXyForIndices(Theta, all_features, matrix, new_cand_theta_idxs(test_idx), new_cand_inst_idxs(test_idx));
% [ypred, ypredvar] = applyModel(model, Xtest, 0, 0, 0);
[rmse, ll, cc, cc_rank] = measures_of_fit(log10(ytest), ypred, ypredvar);