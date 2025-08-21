function algo_config_id = get_algo_config_id(func, x)

K = size(x,2);
switch func.env.algo
    case 'saps'
        table = 'SAPS_CONFIGURATION';
    case 'spear'
%        table = 'SPEAR_CONFIGURATION_1_2';
        table = 'new_SPEAR_CONFIGURATION_1_2';
    case 'cplex'
        table = 'CPLEX_CONFIGURATION_10_1_1';
    case 'cplex12_milp'
        table = 'CPLEX_CONFIGURATION_12_1_MILP';
    case 'cplex12_miqcp'
        table = 'CPLEX_CONFIGURATION_12_1_MIQCP';
    case 'gurobi'
        table = 'GUROBI_CONFIGURATION_2_0_1';
    case 'lpsolve'
        table = 'LPSOLVE_CONFIGURATION_5_5';

    otherwise
        errstr = strcat(['Dont have database entry for algo yet: ', func.env.algo]);
        error(errstr);
end

%=== Backtransform continuous parameters.
x = config_back_transform(x, func);


% %=== Deal with conditional variables: fill irrelevant values with NaN
while 1
    done = 1;
    for i=1:length(func.cond_params_idxs)
        j = func.cond_params_idxs(i);
        if isnan(x(j)), continue, end;
        k = func.parent_param_idxs(i);
        ok_vals = func.ok_parent_value_idxs{i};
        if ~ismember(x(k), ok_vals)
            x(j) = NaN;
            done = 0;
        end
    end
    if done, break, end;
end

%=== Build constraint string.
constraint_string = '';

if strcmp(func.env.algo, 'spear')  && isempty(func.cat)
    constraint_string = 'sp_var_dec_heur is NULL and sp_learned_clause_sort_heur is NULL and sp_orig_clause_sort_heur is NULL and sp_res_order_heur is NULL and sp_clause_del_heur is NULL and sp_phase_dec_heur is NULL and sp_resolution is NULL and ';
end

for i=1:func.dim
    if i>1
        constraint_string = strcat([constraint_string, ' and ']);
    end
    if isnan(x(i))
        value = 'IS NULL';
    else
        if ismember(i, func.cat)
            values = func.all_values{i};
            value = strcat(['= "', values{x(i)}, '"']);
        else
            value = strcat(['= "', num2str(x(i)), '"']);
% the following is now done as part of config_back_transform            
%         else
%             %=== Round continuous parameters to 4 digits -- otherwise DB problems !
%             value = strcat(['= ', num2str(floor(x(i).*10000 + 0.5)/10000)]);
        end
    end
    db_paramname = get_db_paramname(func.param_names{i}, func.env.algo);
    constraint_string = strcat([constraint_string ' ' db_paramname ' ' value]);
end

cmd = strcat(['select ID from ' table ' where ' constraint_string]);
algo_config_id = mysql(cmd);
if length(algo_config_id) > 1 % repeated entries are possible since we can't have an index over ALL parameters (too many)
    algo_config_id = min(algo_config_id);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function paramname = get_db_paramname(paramname, algo)
%=== Replace - by _
dash_idxs = regexp(paramname, '-');
paramname(dash_idxs) = '_';

%=== Fix names for SAPS.
if strcmp(algo, 'saps')
    switch paramname
        case 'alpha'
            paramname = 'ALPHA';
        case 'rho'
            paramname = 'RHO';
        case 'ps'
            paramname = 'SMOOTH_PROB';
        case 'wp'
            paramname = 'RANDOM_WALK_PROB';
        otherwise
            error 'Unkown SAPS parameter.'
    end
end

if strcmp(algo, 'lpsolve')
    if strcmp(paramname, 'Bb')
        paramname = 'Bb_lower';
    elseif strcmp(paramname, 'Bg')
        paramname = 'Bg_lower';
    end
end