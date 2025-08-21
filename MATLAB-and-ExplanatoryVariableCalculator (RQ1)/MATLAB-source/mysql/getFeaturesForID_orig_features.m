function [X, featureNames] = getFeaturesForID(mip, ids)

uniq_ids = my_unique(ids);

global feature_table;
if ~mip
    feature_table = 'FH_THESIS_FEATURES';
%    feature_table = 'NEWFEATURES';
end


for i=1:length(uniq_ids)
    i;
    poss_ids = uniq_ids{i};
    for j=1:length(poss_ids)
        %=== If any possible ID has features computed use those.
        id = poss_ids(j);
        cmd = strcat(['select count(*) from ', feature_table, ' where INSTANCE_ID=', num2str(id)]);
        features_in_db = mysql(cmd);
        if features_in_db
            if mip
                %                 featureNames = {'bg_abs_price_dev'; 'bg_ppg_dev'; 'bg_psqrg_dev'; 'bg_edge_dense'; 'bg_deg_dev'; 'bg_max_deg'; 'bg_min_deg'; 'bp_avg_good_deg'; 'bp_good_deg_dev'; 'bp_max_good_deg'; 'bp_min_good_deg'; 'bp_avg_bid_deg'; 'bp_bid_deg_dev'; 'bp_max_bid_deg'; 'bp_min_bid_deg'; 'bg_deg_q1'; 'bg_deg_q2'; 'bg_deg_q3'; 'lp_avg'; 'lp_l2_avg'; 'lp_linf'; 'lp_objval'; 'cc_bg_radius'; 'cc_bg_ecc_avg'; 'cc_bg_ecc_dev'; 'cc_clustering'; 'cc_path'; 'cc_cp_ratio'; 'cc_clust_dev'};
                %                 cmd = strcat(['select ', featureNames{1}]);
                %                 for j = 2:length(featureNames)
                %                     cmd = strcat([cmd, ', ', featureNames{j}]);
                %                 end
                %                 cmd = strcat([cmd, ' from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]); %not using the following since they're sometimes NULL: bg_avg_deg lp_l1 lp_l2 cc_bg_diameter])
                %                 [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25), tmp(26), tmp(27), tmp(28), tmp(29)] = mysql(cmd);
                featureNames ={'n_constraints'; 'n_variables'; 'nzcnt'; 'probtype'; 'num_q_constr'; 'num_quad'; 'num_qpnz'; 'vcg_var_deg_mean'; 'vcg_var_deg_std'; 'vcg_var_deg_min'; 'vcg_var_deg_max'; 'vcg_con_deg_mean'; 'vcg_con_deg_std'; 'vcg_con_deg_min'; 'vcg_con_deg_max'; 'prices_stddev'; 'price_per_number_stddev'; 'price_per_sqrt_number_stddev'; 'perc_cont_var'; 'support_size_mean'; 'support_size_std'; 'a_ij_magnitude_mean'; 'a_ij_magnitude_std'; 'a_ij_varcoef_avg'; 'a_ij_varcoef_std'; 'perc_unbounded_discrete'; 'edge_density'; 'vg_deg_std'; 'vg_deg_max'; 'vg_deg_min'; 'vg_deg_q0_25'; 'vg_deg_med'; 'vg_deg_q0_75'; 'clust_coef'; 'deviation'; 'lp_avg'; 'lp_l2_avg'; 'lp_linf'; 'lp_objval'};
                cmd = strcat(['select ', featureNames{1}]);
                for j = 2:length(featureNames)
                    cmd = strcat([cmd, ', ', featureNames{j}]);
                end
                cmd = strcat([cmd, ' from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);
                [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25), tmp(26), tmp(27), tmp(28), tmp(29), tmp(30), tmp(31), tmp(32), tmp(33), tmp(34), tmp(35), tmp(36), tmp(37), tmp(38), tmp(39)] = mysql(cmd);
            else
                if 0

                    %            cmd = strcat(['select saps_BestSolution_CoeffVariance*saps_BestSolution_CoeffVariance, saps_BestSolution_CoeffVariance*saps_AvgImproveToBest_Mean, saps_BestStep_CoeffVariance*saps_FirstLMRatio_Mean, gsat_BestSolution_CoeffVariance*lobjois_mean_depth_over_vars, saps_AvgImproveToBest_CoeffVariance, saps_BestCV_Mean*gsat_BestStep_Mean, VG_mean*gsat_BestStep_Mean, saps_AvgImproveToBest_CoeffVariance*gsat_BestSolution_Mean, vars_clauses_ratio*lobjois_mean_depth_over_vars, VG_mean*gsat_BestStep_Mean, TRINARY_PLUS*lobjois_log_num_nodes_over_vars from ', feature_table, ' where INSTANCE_ID=', num2str(id)]);
                    %    cmd = strcat('select saps_BestSolution_CoeffVariance*saps_BestSolution_CoeffVariance, saps_BestSolution_CoeffVariance*saps_AvgImproveToBest_Mean, saps_BestStep_CoeffVariance*saps_FirstLMRatio_Mean, gsat_BestSolution_CoeffVariance*lobjois_mean_depth_over_vars, saps_AvgImproveToBest_CoeffVariance, saps_BestCV_Mean*gsat_BestStep_Mean, VG_mean*gsat_BestStep_Mean, saps_AvgImproveToBest_CoeffVariance*gsat_BestSolution_Mean, vars_clauses_ratio*lobjois_mean_depth_over_vars, VG_mean*gsat_BestStep_Mean, TRINARY_PLUS*lobjois_log_num_nodes_over_vars from FH_FEATURES_1sLS where INSTANCE_ID=', num2str(id));
                    %            [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11)] = mysql(cmd);

                    % "Original features" from CP-06
                    % cmd = strcat(['select nvars, nclauses, vars_clauses_ratio, VCG_VAR_min, VCG_VAR_max, VCG_VAR_entropy, VCG_CLAUSE_mean, VCG_CLAUSE_coeff_variation, VCG_CLAUSE_max, VCG_CLAUSE_entropy, POSNEG_RATIO_CLAUSE_mean, POSNEG_RATIO_VAR_stdev, POSNEG_RATIO_VAR_max, BINARY_PLUS, TRINARY_PLUS, HORNY_VAR_coeff_variation, HORNY_VAR_min, HORNY_VAR_max, HORNY_VAR_entropy, horn_clauses_fraction, VG_mean, VG_coeff_variation, VG_min, vars_reduced_depth_4, vars_reduced_depth_16, vars_reduced_depth_64, saps_BestSolution_Mean, saps_BestSolution_CoeffVariance, saps_BestStep_Mean, saps_BestStep_CoeffVariance, saps_AvgImproveToBest_Mean, saps_AvgImproveToBest_CoeffVariance, saps_FirstLMRatio_Mean, saps_FirstLMRatio_CoeffVariance, saps_BestCV_Mean, gsat_BestSolution_Mean, gsat_BestSolution_CoeffVariance, gsat_BestStep_Mean, gsat_BestStep_CoeffVariance, gsat_AvgImproveToBest_CoeffVariance, gsat_BestCV_Mean, lobjois_mean_depth_over_vars, lobjois_log_num_nodes_over_vars from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);

                    %taken 3 out b/c very costly to compute: VG_mean, VG_coeff_variation,
                    %VG_min,
                    featureNames = {'nvars'; 'nclauses'; 'vars_clauses_ratio'; 'VCG_VAR_min'; 'VCG_VAR_max'; 'VCG_VAR_entropy'; 'VCG_CLAUSE_mean'; 'VCG_CLAUSE_coeff_variation'; 'VCG_CLAUSE_max'; 'VCG_CLAUSE_entropy'; 'POSNEG_RATIO_CLAUSE_mean'; 'POSNEG_RATIO_VAR_stdev'; 'POSNEG_RATIO_VAR_max'; 'BINARY_PLUS'; 'TRINARY_PLUS'; 'HORNY_VAR_coeff_variation'; 'HORNY_VAR_min'; 'HORNY_VAR_max'; 'HORNY_VAR_entropy'; 'horn_clauses_fraction'; 'vars_reduced_depth_4'; 'vars_reduced_depth_16'; 'vars_reduced_depth_64'; 'lobjois_mean_depth_over_vars'; 'lobjois_log_num_nodes_over_vars'};
                    cmd = strcat(['select ', featureNames{1}]);
                    for j = 2:length(featureNames)
                        cmd = strcat([cmd, ', ', featureNames{j}]);
                    end
                    cmd = strcat([cmd, ' from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);
                    %    [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25), tmp(26), tmp(27), tmp(28)] = mysql(cmd);
                    [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25)] = mysql(cmd);

                    lsFeatureNames = {'saps_BestSolution_Mean'; 'saps_BestSolution_CoeffVariance'; 'saps_BestStep_Mean'; 'saps_BestStep_CoeffVariance'; 'saps_AvgImproveToBest_Mean'; 'saps_AvgImproveToBest_CoeffVariance'; 'saps_FirstLMRatio_Mean'; 'saps_FirstLMRatio_CoeffVariance'; 'saps_BestCV_Mean'; 'gsat_BestSolution_Mean'; 'gsat_BestSolution_CoeffVariance'; 'gsat_BestStep_Mean'; 'gsat_BestStep_CoeffVariance'; 'gsat_AvgImproveToBest_CoeffVariance'; 'gsat_BestCV_Mean'};
                    cmd = strcat(['select ', lsFeatureNames{1}]);
                    for j = 2:length(lsFeatureNames)
                        cmd = strcat([cmd, ', ', lsFeatureNames{j}]);
                    end
                    %    cmd = strcat([cmd, ' from lsFEATURES where INSTANCE_ID = ', num2str(id)]);
                    cmd = strcat([cmd, ' from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);
                    [tmp(26), tmp(27), tmp(28), tmp(29), tmp(30), tmp(31), tmp(32), tmp(33), tmp(34), tmp(35), tmp(36), tmp(37), tmp(38), tmp(39), tmp(40)] = mysql(cmd);
                    %    [tmp(29), tmp(30), tmp(31), tmp(32), tmp(33), tmp(34), tmp(35), tmp(36), tmp(37), tmp(38), tmp(39), tmp(40), tmp(41), tmp(42), tmp(43)] = mysql(cmd);

                    featureNames = [featureNames; lsFeatureNames];
                    %=== "New" set of 48 features we use for SATzilla. Lin says lots of the
                    %=== previously used local search features were unstable.
                    %=== Minus the 4 VG features as they are sometimes super expensive:
                    %=== VG_mean, VG_coeff_variation, VG_min, VG_max. Without them nothing takes longer than 5 seconds!

                    %                 cmd = strcat(['select nvars, nclauses, vars_clauses_ratio, VCG_VAR_mean, VCG_VAR_coeff_variation, VCG_VAR_min, VCG_VAR_max, VCG_VAR_entropy, VCG_CLAUSE_mean, VCG_CLAUSE_coeff_variation, VCG_CLAUSE_min, VCG_CLAUSE_max, VCG_CLAUSE_entropy, POSNEG_RATIO_CLAUSE_mean, POSNEG_RATIO_CLAUSE_coeff_variation, POSNEG_RATIO_CLAUSE_entropy, POSNEG_RATIO_VAR_mean, POSNEG_RATIO_VAR_stdev, POSNEG_RATIO_VAR_min, POSNEG_RATIO_VAR_max, POSNEG_RATIO_VAR_entropy, BINARY_PLUS, TRINARY_PLUS, HORNY_VAR_mean, HORNY_VAR_coeff_variation, HORNY_VAR_min, HORNY_VAR_max, HORNY_VAR_entropy, horn_clauses_fraction, vars_reduced_depth_1, vars_reduced_depth_4, vars_reduced_depth_16, vars_reduced_depth_64, vars_reduced_depth_256, lobjois_mean_depth_over_vars, lobjois_log_num_nodes_over_vars from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);
                    %                 [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25), tmp(26), tmp(27), tmp(28), tmp(29), tmp(30), tmp(31), tmp(32), tmp(33), tmp(34), tmp(35), tmp(36)] = mysql(cmd);
                    %
                    % %                cmd = strcat(['select saps_BestStep_Mean, saps_BestStep_Median, saps_BestStep_Q10, saps_BestStep_Q90, saps_AvgImproveToBest_Mean, saps_FirstLMRatio_Mean, saps_BestCV_Mean, gsat_FirstLMRatio_Mean from lsFEATURES where INSTANCE_ID = ', num2str(id)]);
                    %                 cmd = strcat(['select saps_BestStep_Mean, saps_BestStep_Median, saps_BestStep_Q10, saps_BestStep_Q90, saps_AvgImproveToBest_Mean, saps_FirstLMRatio_Mean, saps_BestCV_Mean, gsat_FirstLMRatio_Mean from FEATURES where INSTANCE_ID = ', num2str(id)]);
                    %                 [tmp(37), tmp(38), tmp(39), tmp(40), tmp(41), tmp(42), tmp(43), tmp(44)] = mysql(cmd);
                else

                    featureNames = {'nvars', 'nclauses', 'vars_clauses_ratio', 'VCG_VAR_mean', 'VCG_VAR_coeff_variation', 'VCG_VAR_min', 'VCG_VAR_max', 'VCG_VAR_entropy', 'VCG_CLAUSE_mean', 'VCG_CLAUSE_coeff_variation', 'VCG_CLAUSE_min', 'VCG_CLAUSE_max', 'VCG_CLAUSE_entropy', 'POSNEG_RATIO_CLAUSE_mean', 'POSNEG_RATIO_CLAUSE_coeff_variation', 'POSNEG_RATIO_CLAUSE_min', 'POSNEG_RATIO_CLAUSE_max', 'POSNEG_RATIO_CLAUSE_entropy', 'POSNEG_RATIO_VAR_mean', 'POSNEG_RATIO_VAR_stdev', 'POSNEG_RATIO_VAR_min', 'POSNEG_RATIO_VAR_max', 'POSNEG_RATIO_VAR_entropy', 'UNARY', 'BINARY_PLUS', 'TRINARY_PLUS', 'HORNY_VAR_mean', 'HORNY_VAR_coeff_variation', 'HORNY_VAR_min', 'HORNY_VAR_max', 'HORNY_VAR_entropy', 'horn_clauses_fraction', 'VG_mean', 'VG_coeff_variation', 'VG_min', 'VG_max', 'CG_mean', 'CG_coeff_variation', 'CG_min', 'CG_max', 'CG_entropy', 'cluster_coeff_mean', 'cluster_coeff_coeff_variation', 'cluster_coeff_min', 'cluster_coeff_max', 'cluster_coeff_entropy', 'vars_reduced_depth_1', 'vars_reduced_depth_4', 'vars_reduced_depth_16', 'vars_reduced_depth_64', 'vars_reduced_depth_256', 'saps_BestSolution_Mean', 'saps_BestSolution_CoeffVariance', 'saps_BestStep_Mean', 'saps_BestStep_CoeffVariance', 'saps_BestStep_Median', 'saps_BestStep_Q10', 'saps_BestStep_Q90', 'saps_AvgImproveToBest_Mean', 'saps_AvgImproveToBest_CoeffVariance', 'saps_FirstLMRatio_Mean', 'saps_FirstLMRatio_CoeffVariance', 'saps_BestCV_Mean', 'gsat_BestSolution_Mean', 'gsat_BestSolution_CoeffVariance', 'gsat_BestStep_Mean', 'gsat_BestStep_CoeffVariance', 'gsat_BestStep_Median', 'gsat_BestStep_Q10', 'gsat_BestStep_Q90', 'gsat_AvgImproveToBest_Mean', 'gsat_AvgImproveToBest_CoeffVariance', 'gsat_FirstLMRatio_Mean', 'gsat_FirstLMRatio_CoeffVariance', 'gsat_BestCV_Mean', 'lobjois_mean_depth_over_vars', 'lobjois_log_num_nodes_over_vars'};
                    cmd = strcat(['select ', featureNames{1}]);
                    for j = 2:length(featureNames)
                        cmd = strcat([cmd, ', ', featureNames{j}]);
                    end
                    cmd = strcat([cmd, ' from ', feature_table, ' where INSTANCE_ID = ', num2str(id)]);
                    [tmp(1), tmp(2), tmp(3), tmp(4), tmp(5), tmp(6), tmp(7), tmp(8), tmp(9), tmp(10), tmp(11), tmp(12), tmp(13), tmp(14), tmp(15), tmp(16), tmp(17), tmp(18), tmp(19), tmp(20), tmp(21), tmp(22), tmp(23), tmp(24), tmp(25), tmp(26), tmp(27), tmp(28), tmp(29), tmp(30), tmp(31), tmp(32), tmp(33), tmp(34), tmp(35), tmp(36), tmp(37), tmp(38), tmp(39), tmp(40), tmp(41), tmp(42), tmp(43), tmp(44), tmp(45), tmp(46), tmp(47), tmp(48), tmp(49), tmp(50), tmp(51), tmp(52), tmp(53), tmp(54), tmp(55), tmp(56), tmp(57), tmp(58), tmp(59), tmp(60), tmp(61), tmp(62), tmp(63), tmp(64), tmp(65), tmp(66), tmp(67), tmp(68), tmp(69), tmp(70), tmp(71), tmp(72), tmp(73), tmp(74), tmp(75), tmp(76), tmp(77)] = mysql(cmd);

                end
            end
            uniqX(i,:) = tmp;
            break
        else
            str = strcat(['No features exist for poss_ids ', num2str(poss_ids)]);
            error(str)
        end
    end
end

X = [];
for i = 1:length(ids)
    for j=1:length(uniq_ids)
        %=== It's enough if the first array element is the same.
        if ids{i}(1) == uniq_ids{j}(1)
            X(i,:) = uniqX(j,:);
            break;
        end
    end
end


function uniq_cell = my_unique(cell_with_arrays)
uniq_cell = {};
for i=1:length(cell_with_arrays)
    contained = 0;
    for j=1:length(uniq_cell)
        %=== It's enough if the first array element is the same.
        if cell_with_arrays{i}(1) == uniq_cell{j}(1)
            contained = 1;
            break;
        end
    end
    if ~contained
        uniq_cell{end+1} = cell_with_arrays{i};
    end
end