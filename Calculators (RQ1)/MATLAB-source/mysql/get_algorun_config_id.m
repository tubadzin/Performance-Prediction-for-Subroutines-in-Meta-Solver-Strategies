function algorun_config_id = get_algorun_config_id(algo, algoconfig_id, inst_id)
cmd = strcat('select ID from AL_ALGORUN_CONFIG where ALGORITHM = "', algo, '" and INSTANCE_ID = ', num2str(inst_id), ' and ALGORITHM_CONFIG_ID = ', num2str(algoconfig_id));
algorun_config_id = mysql(cmd);

if isempty(algorun_config_id)
    cmd = strcat('select ID from MORE_AL_ALGORUN_CONFIG where ALGORITHM = "', algo, '" and INSTANCE_ID = ', num2str(inst_id), ' and ALGORITHM_CONFIG_ID = ', num2str(algoconfig_id));
    algorun_config_id = mysql(cmd);
end