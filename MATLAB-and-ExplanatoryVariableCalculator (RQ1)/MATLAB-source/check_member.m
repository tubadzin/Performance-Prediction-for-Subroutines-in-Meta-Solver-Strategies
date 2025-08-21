function result = check_member(value, cell_array)
result = (value <= length(cell_array));
% for i=1:length(cell_array)
%     if strcmp(cell_array{i}, num2str(value))
%         result = true;
%         return;
%     end
% end
% result = false;