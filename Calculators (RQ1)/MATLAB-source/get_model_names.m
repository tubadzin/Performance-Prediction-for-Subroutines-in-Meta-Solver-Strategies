function model_names = get_model_names(options_vec, fixname)
for i=1:length(options_vec)
    name = options_vec{i}.unique_model_name;
%     name = options_vec{i}.modelType;
%     switch name
%         case 'rf' 
%             if isfield(options_vec{i}, 'opt') && options_vec{i}.opt
%                 name = 'RF cv';
%             else
%                 name = 'RF';
%             end
%         case 'regression-tree'
%             name = 'RT';
%         case 'GPML'
%             name = 'PP';
%         case 'LR'
%             name = 'RR';
%     end
    if fixname
        model_names{i} = fix_modelname(name);
    else
        model_names{i} = name;
    end
end