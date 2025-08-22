function algo_config_id = get_saps_config_id(alpha, rho, ps, wp)
error 'just use the general one'

assert(length(alpha)==1, 'DB functions not vectorized.');

cmd = strcat('select ID from SAPS_CONFIGURATION where ALPHA = ', num2str(alpha), ' and RHO = ', num2str(rho), ' and SMOOTH_PROB = ', num2str(ps), ' and RANDOM_WALK_PROB = ', num2str(wp));
algo_config_id = mysql(cmd);