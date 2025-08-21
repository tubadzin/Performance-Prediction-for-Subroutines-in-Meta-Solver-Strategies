function samples = fit_dist_and_sample(y_pred_all, cens_pred_all, weights_all, numSamples, lowerBoundForSamples, valueForAllCens)

%=== Prediction with a KM estimator, taking its mean as the area under the curve.
uncens_idx = find(cens_pred_all==0);
cens_idx = find(cens_pred_all==1);
y = y_pred_all(uncens_idx);
c = y_pred_all(cens_idx);

%=== Fit KM to it, compute its mean as area under the curve.
[km_mean, km_median, t, F] = km_stats_w(y, c, valueForAllCens, weights_all(uncens_idx), weights_all(cens_idx));

%=== Get quantile of lowerBoundForSamples, and only sample above that
%quantile.

% u=0:1/(numSamples+1):1;
% u = u(2:end-1);
% if isempty(y)
if isempty(y) || max(y) < lowerBoundForSamples-1e-10
%    warning 'Only censored data points - using upper bound as fit';
    %error 'Cannot fit distribution to only censored data points!'
    samples(1:numSamples) = valueForAllCens;
else
    %=== Get quantile of lowerBoundForSamples, and only sample above that
    %quantile.
    if lowerBoundForSamples+1e-10 < t(1)
        lb_quantile = 0;
    end
    for i=2:length(t)
        if t(i) >= lowerBoundForSamples+1e-10
            lb_quantile = 1-F(i-1);
            break
        else
            lb_quantile = 1-F(i);
        end
    end
    u = lb_quantile:(1-lb_quantile)/(numSamples+1):1;
    u = u(2:end-1);
    
    km_max = 1-F(end);
    u_km = u(find(u<=km_max));
    u_par = u(find(u>km_max));

    samples = zeros(1,numSamples);
    for j=1:length(u_km)
        idx = find(F <= 1-u_km(j) + 1e-6);
        samples(j) = t(idx(1));
    end

    if ~isempty(u_par)
%             %=== Fit extreme value distribution to it.
%         %         [parmhat, parmci] = evfit([1e-5,mean(10.^y_train),y_pred_all],0.05,[0,0,cens_pred_all],[1e-5, 1e-5, ones(1,length(y_pred_all))]);
%             [parmhat, parmci] = evfit(y_pred_all,0.05,cens_pred_all,weights_all);
%             if any(isnan(parmhat))
%                 parmhat = parmhat;
%                 number_of_uncensored = length(find(cens_pred_all==0));
%                 number_of_censored = length(find(cens_pred_all==1));
% 
%                 %errstr = strcat(['parmhat has NaN entries -- cannot fit Weibull to ', num2str(number_of_uncensored), ' uncensored points and ', num2str(number_of_censored), ' censored ones.']);
%                 %error(errstr);
%                 fprintf('WARNING: cannot fit extreme value distribution to ', num2str(number_of_uncensored), ' uncensored points and ', num2str(number_of_censored), ' censored ones. Using default delta distribution cutoff+1.');
%                 [parmhat, parmci] = evfit([upper_mean]); 
%                 sample(length(u_km)+1:length(u)) = upper_mean;
%             else
%                 sample(length(u_km)+1:length(u)) = evinv(u_par, parmhat(1), parmhat(2));
%                 %sample_pred(i, :) = evrnd(parmhat(1), parmhat(2), numSamples,1);
% 
%                 if any(isnan(sample(length(u_km)+1:length(u))) ) % if problems with stratified sampling, just sample directly.
%                     sample(length(u_km)+1:length(u)) = evrnd(parmhat(1), parmhat(2), numSamples,1);
%                 end
%             end

    %=== Fit Weibull to it, get samples from that.
        [parmhat, parmci] = wblfit(y_pred_all,0.05,cens_pred_all,weights_all);
        if any(isnan(parmhat))
            parmhat = parmhat
            number_of_uncensored = length(find(cens_pred_all==0));
            number_of_censored = length(find(cens_pred_all==1));

            %errstr = strcat(['parmhat has NaN entries -- cannot fit Weibull to ', num2str(number_of_uncensored), ' uncensored points and ', num2str(number_of_censored), ' censored ones.']);
            %error(errstr);
            fprintf('WARNING: cannot fit Weibull to ', num2str(number_of_uncensored), ' uncensored points and ', num2str(number_of_censored), ' censored ones. Using default delta distribution at valueForAllCens.');
%             [parmhat, parmci] = evfit([upper_mean]); 
            samples(length(u_km)+1:length(u)) = valueForAllCens;
        else
            samples(length(u_km)+1:length(u)) = wblinv(u_par, parmhat(1), parmhat(2));
            %sample_pred(i, :) = wblrnd(parmhat(1), parmhat(2), numSamples,1);

            if any(isnan(samples(length(u_km)+1:length(u))) ) % if problems with stratified sampling, just sample directly.
                warning 'stratified sampling did not work - sampling directly.';
                samples(length(u_km)+1:length(u)) = wblrnd(parmhat(1), parmhat(2), numSamples,1);
            end
        end
    end
    samples = samples(randperm(numSamples)); % randperm such that we don't have the same u for each data point in a tree
end
assert(all(~isinf(samples)));