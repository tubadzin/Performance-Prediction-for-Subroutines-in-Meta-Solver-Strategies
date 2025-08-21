function [theta_idxs, inst_idxs] = get_idx_pairs(theta_idx_candidates, inst_idx_candidates)
theta_idxs = -ones(length(theta_idx_candidates)*length(inst_idx_candidates), 1);
inst_idxs = -ones(length(theta_idx_candidates)*length(inst_idx_candidates), 1);
count = 1;
for i=theta_idx_candidates
    for j=inst_idx_candidates
        theta_idxs(count) = i;
        inst_idxs(count) = j;
        count = count+1;
    end
end
assert(all(theta_idxs>0));
assert(all(inst_idxs>0));