function header_output(modelNames, statsNames, extra1, extra2)
if nargin < 3
    extra1 = '';
    extra2 = '';
end

N = length(statsNames);

str = strcat('Domain & ', extra1);
for i=1:N
    if i < N
        str = strcat([str, ' \\multicolumn{', num2str(length(modelNames)), '}{|c|}{', statsNames{i}, '}']);
        str = strcat([str, ' & ']);
    else
        str = strcat([str, ' \\multicolumn{', num2str(length(modelNames)), '}{|c}{', statsNames{i}, '}']);
        str = strcat(str, '\\\\\n');
    end
end

str = strcat(str, ' & ', extra2);
for i=1:N
    if i < N
        for j=1:length(modelNames)
            str = strcat([str, fix_modelname(modelNames{j}), ' & ']);
        end
    else
        for j=1:length(modelNames)
            str = strcat([str, fix_modelname(modelNames{j})]);
            if j < length(modelNames)
                str = strcat([str, ' & ']);
            else
                str = strcat(str, '\\\\\n');        
            end
        end
    end
end

fprintf(str);