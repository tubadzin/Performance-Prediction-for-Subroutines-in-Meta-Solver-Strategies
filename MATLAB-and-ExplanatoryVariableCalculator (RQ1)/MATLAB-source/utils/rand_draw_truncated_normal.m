function samples = rand_draw_truncated_normal(a, b, mu, sigma, sampleSize, u)
% Generates samples from the truncated Normal distribution restricted to
% the interval [a,b]. If u in [0,1] is given, then this returns the 100*u-th
% percentile.

if nargin < 6
    u = rand(sampleSize);
end

PHIl = normcdf((a-mu)/sigma);
PHIr = normcdf((b-mu)/sigma);

samples = mu + sigma*( sqrt(2)*erfinv(2*(PHIl+(PHIr-PHIl)*u)-1) );                    