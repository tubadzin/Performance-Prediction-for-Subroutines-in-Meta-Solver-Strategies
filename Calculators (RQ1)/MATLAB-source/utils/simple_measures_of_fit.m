function [rmse_w, ll_w, cc_w, cc_rank] = simple_measures_of_fit(logModel, y, y_pred, y_predvar, cens, weights)
%=== Compute RMSE, CC, and LL; the first two only on noncensored data.
%=== Also compute weighted versions.
if nargin < 6
    weights = ones(length(y),1);
    if nargin < 5
        cens = zeros(length(y),1);
        if nargin < 4
            y_predvar = zeros(length(y),1);
        end
    end
end

if logModel
    y = log10(y);
end
[rmse_w, ll_w, cc_w, cc_rank] = measures_of_fit(y, y_pred, y_predvar, cens, weights);