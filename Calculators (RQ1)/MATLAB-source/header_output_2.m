function header_output_2(modelNames, statsNames)

N = length(modelNames);

str = 'Domain & ';
for i=1:N
    if i < N
        str = strcat([str, ' \\multicolumn{', num2str(length(statsNames)), '}{|c|}{', modelNames{i}, '}']);
        str = strcat([str, ' & ']);
    else
        str = strcat([str, ' \\multicolumn{', num2str(length(statsNames)), '}{|c}{', modelNames{i}, '}']);
        str = strcat(str, '\\\\\n');
    end
end

for i=1:N
    for j=1:length(statsNames)
        str = strcat([str, ' & ', statsNames{j}]);
    end
end
str = strcat(str, '\\\\\n');

fprintf(str);