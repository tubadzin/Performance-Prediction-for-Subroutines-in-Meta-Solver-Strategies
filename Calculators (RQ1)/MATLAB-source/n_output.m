function n_output(values,small_good)

% Bring cell into array format.
if iscell(values)
    values_as_array = [];
    for i=1:length(values)
        values_as_array(i) = values{i};
    end
    values = values_as_array;
end

N = length(values);
values = two_dec(values);
eps = 0.0001;

best_ind = zeros(N,1);
if small_good
    best = min(values);
    for i=1:N
        if values(i) <= best + eps
            best_ind(i) = 1;
        end
    end
else
    best = max(values);
    for i=1:N
        if values(i) >= best - eps
            best_ind(i) = 1;
        end
    end
end

str = '';
for i=1:N
    if best_ind(i)
        str = strcat([str, ' \\textbf{', to_s(values(i)), '}']);
    else
        str = strcat([str, ' ', to_s(values(i))]);
    end
     if i < N
        str = strcat([str, ' & ']);
     end
end
% str = strcat(str, '\n');
fprintf(str);

function x=two_dec(x)
x = ceil(x*100-0.5)/100;

function x = to_s(x)
% x = sprintf( '%4.2f' , x );
%x = sprintf( '%0.2g' , x );
x = sprintf( '%0.2g' , x );
