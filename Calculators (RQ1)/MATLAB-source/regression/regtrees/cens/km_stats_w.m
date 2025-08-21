%function [km_mean, km_median, t, F] = km_stats(y, c)
function [km_mean, km_median, t, F] = km_stats_w(y, c, upper_bound, weights_y, weights_c)
% Kaplan-Meier estimator. If the highest value is censored K-L is
% undefined for values above that.
% For computing the mean we need an estimate, though. Fortunately, we
% have an upper bound, and we assume that there are no deaths until
% that upper bound.

% assert(all(y>=0));
% assert(all(c>=0));

N = length(y);
M = length(c);
if N+M==0
    error 'Kaplan-Meier estimator undefined for empty population.'
end

[c, idx] = sort(c);
[y, idx2] = sort(y);
weights_c = weights_c(idx);
weights_y = weights_y(idx2);

c_idx = 1;
prod=1;
km_mean = 0;
t=[];
F=[];
median_set = 0;
N_i = sum(weights_y) + sum(weights_c);

if N>0
    rows = find(y(1:end-1)+eps < y(2:end));
    rows(end+1) = length(y); % y is now guaranteed to have at last one element.
    for i=1:length(rows)
        t(i) = y(rows(i));
        d_i = 0;
        while c_idx <= length(c) && c(c_idx) < t(i)-1e-6
            N_i = N_i-weights_c(c_idx);
            c_idx = c_idx + 1;
        end
        if i==1
            d_i = sum(weights_y(1:rows(i)));
            km_mean = km_mean + prod * y(rows(i));
        else
            d_i = sum(weights_y(rows(i-1)+1:rows(i)));
            km_mean = km_mean + prod * (y(rows(i)) - y(rows(i-1)));
        end
        prod = prod * (N_i-d_i)/(N_i+0.0);
        N_i = N_i-d_i;
        F(i) = prod;
        if ~median_set && prod < 0.5 + eps
            if prod < 0.5 - eps
                km_median = t(i);
            else
                if i<length(rows)
                    km_median = mean([t(i), y(rows(i+1))]);
                else
                    km_median = mean([t(i), upper_bound]);
                end
            end
            median_set = 1;
       end
    end
    %=== Deal with remaining censored values (if the highest value is
    %=== uncensored, then prod=0, so nothing happens then)
%    km_mean = km_mean + prod * (upper_bound-y(end));
    if prod > 1e-6
        km_mean = km_mean + prod * (upper_bound-y(end))/2; % uniform distribution between highest and upper bound, mean in the middle.
    end
else
    km_mean = upper_bound;
end

if ~median_set
    km_median = upper_bound;
end