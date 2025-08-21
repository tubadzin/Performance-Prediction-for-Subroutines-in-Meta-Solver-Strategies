function inst_names = get_inst_name(inst_ids)
inst_names = cell(length(inst_ids),1);
for i=1:length(inst_ids)
    inst_id = inst_ids(i);
    cmd = strcat('select concat(RELATIVE_PATH, "/", FILENAME) from INSTANCE where INSTANCE_ID = ', num2str(inst_id));
    inst_names{i} = mysql(cmd);
end