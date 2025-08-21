clear
clear mex
clear functions
rehash path
rehash toolbox
rehash toolboxcache
which -all fh_random_regtreefit_big_leaves_twofeaturetypes_dist_partition


project_root = '/Users/piotrmalkowski/Desktop/Bachelor/Hutter et al./source';
addpath(genpath(project_root))   

ghostDir   = '/opt/homebrew/bin';

setenv('PATH',[getenv('PATH') ':' ghostDir]);   % Matlab now finds 'gs'
setenv('GHOSTSCRIPT',[ghostDir '/gs']);         % export_fig looks here first

!gs -v


% epm_experiments(1,1)





%cd('/Users/piotrmalkowski/Desktop/Bachelor/Hutter et al./source/regression/regtrees')

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%fh_random_regtreefit_big_leaves_twofeaturetypes_dist_partition.c

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%fwd_big_var.c

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%fwd_big.c

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%compute_from_leaves_part.c

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%collect_big_leaves_theta_pis_distinleaf_nomissing.c

%mex -v -compatibleArrayDims ...
%CFLAGS="\$CFLAGS -O3 -std=c99" ... 
%cens_fh_random_regtreefit_big_leaves_twofeaturetypes_dist_part.c