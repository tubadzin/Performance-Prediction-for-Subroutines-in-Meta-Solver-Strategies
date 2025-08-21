function inst_ids = get_inst_id(inst_names)
% For each filename return all matching instances -- necessary because of
% stupid issues in adding instances to the database.

uniq_inst_names = unique(inst_names);

%=== Get the ids for all unique filenames.
uniq_inst_ids = cell(length(uniq_inst_names),1);
for i=1:length(uniq_inst_names)
    poss_ids = [];
    inst_name = uniq_inst_names{i};
    idx = strfind(inst_name, '/');
    path = inst_name(1:idx(end)-1);
    filename = inst_name(idx(end)+1:end);

    cmd = strcat(['select INSTANCE_ID from INSTANCE where RELATIVE_PATH = "' path '" and FILENAME = "' filename '"']);
    res = mysql(cmd);
    if ~isempty(res)
        poss_ids(end+1) =res;
    end

    orig_filename = filename;
    filename = strcat([filename, '\n']);
    cmd = strcat(['select INSTANCE_ID from INSTANCE where RELATIVE_PATH = "' path '" and FILENAME = "' filename '"']);
    res = mysql(cmd);
    if ~isempty(res)
        poss_ids(end+1) =res;
    end
    
    path = strcat([path, '/']);
    cmd = strcat(['select INSTANCE_ID from INSTANCE where RELATIVE_PATH = "' path '" and FILENAME = "' filename '"']);
    res = mysql(cmd);
    if ~isempty(res)
        poss_ids(end+1) =res;
    end
    
    filename = orig_filename;
    cmd = strcat(['select INSTANCE_ID from INSTANCE where RELATIVE_PATH = "' path '" and FILENAME = "' filename '"']);
    res = mysql(cmd);
    if ~isempty(res)
        poss_ids(end+1) =res;
    end
    
    uniq_inst_ids{i} = poss_ids;
end

%=== Return the ids for the filenames in the original order, including repetitions
inst_ids = cell(length(inst_names), 1);
for i=1:length(inst_names)
    for j=1:length(uniq_inst_names)
        if strcmp(inst_names{i}, uniq_inst_names{j})
            inst_ids{i} = uniq_inst_ids{j};
            break;
        end
    end
end