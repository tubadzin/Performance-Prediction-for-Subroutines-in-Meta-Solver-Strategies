function [runtime, censored, runlength, best_sol, solved] = get_runtime(algorun_config_id, seed, cutoff_time, dorunlength)
% OLD DB: cmd = strcat(['select MEASURED_RUNTIME from AL_ALGORUN where ALGORUN_CONFIG_ID = ' num2str(algorun_config_id) ' and SEED = ' num2str(seed) ' and CUTOFF_TIME = ' num2str(cutoff_time)]);

if dorunlength
    cutoff_length = cutoff_time; 
    cutoff_time = 5;% HACK since CUTOFF_LENGTH not saved in DB -> in this time we have to do more steps than CUTOFF_LENGTH.
else
    cutoff_length = 1000000;
end
cmd = strcat(['select SOLVED, BEST_SOL, MEASURED_RUNTIME, MEASURED_RUNLENGTH, CUTOFF_TIME from AL_ALGORUN where ALGORUN_CONFIG_ID = ' num2str(algorun_config_id) ' and SEED = ' num2str(seed)]);
[solveds, best_sols, runtimes, runlengths, cutoff_times] = mysql(cmd);
    
cmd = strcat(['select SOLVED, BEST_SOL, MEASURED_RUNTIME, MEASURED_RUNLENGTH, CUTOFF_TIME from MORE_AL_ALGORUN where ALGORUN_CONFIG_ID = ' num2str(algorun_config_id) ' and SEED = ' num2str(seed)]);
[more_solveds, more_best_sols, more_runtimes, more_runlengths, more_cutoff_times] = mysql(cmd);
solveds = [solveds; more_solveds];
best_sols = [best_sols; more_best_sols];
runtimes = [runtimes; more_runtimes];
runlengths = [runlengths; more_runlengths];
cutoff_times = [cutoff_times; more_cutoff_times];

for i=1:length(solveds)
    
    if strcmp(solveds{i}, 'TIMEOUT')
        solved = 0;
    elseif strcmp(solveds{i}, 'SAT')
        solved = 1;                    
    elseif strcmp(solveds{i}, 'UNSAT')
        solved = 2;
    elseif strcmp(solveds{i}, 'WRONG ANSWER') || strcmp(solveds{i}, 'WRONG')
        solved = 3;
%         runtime = inf;
%         runlength = inf;
%         best_sol = inf;
    else
        error 'unknown solution status'
    end

    if solved == 3
        warnstring = strcat(['Wrong answer for algorun_config_id ', num2str(algorun_config_id), ' and seed ', num2str(seed)]);
        warning(warnstring);
        bout(warnstring);
        
        censored = 0;
        runtime = inf;
        runlength = inf;
        solved = 3;
        best_sol = inf;
%         runtime = runtimes(i);
%         runlength = runlengths(i);
%         censored = 1;
%         best_sol = best_sols(i);

        return;
        
    elseif solved == 0
        if (~dorunlength && runtimes(i) < cutoff_times(i))
            %=== Spear cannot deal with non-integer captimes => problems in
            %the database that we now have to fix...
            runtimes
            cutoff_times
            solveds
            warning('(was warning before: timeout but measured_runtime < cutoff_time)')

            censored = 1;
            runtime = cutoff_times(i) + 0.001;
            runlength = cutoff_length;
            solved = 0;
            best_sol = best_sols(i);
            
            return;
        end
        if (dorunlength && runlengths(i) >= cutoff_length) || (~dorunlength && cutoff_times(i) >= cutoff_time - 0.0001 && cutoff_times(i) <= cutoff_time + 0.0001)
            %=== A run with exactly this cutoff timed out => timeout
            censored = 1;
            runtime = cutoff_time;
            runlength = cutoff_length;
            solved = 0;
            best_sol = best_sols(i);
            return
        end
% This case can only return the solution quality with a larger kappa; that
% may lead to problems when we optimize solution quality for different
% captimes. => rather redo the run in this case than using the cached run.
%         if (dorunlength && runlengths(i) >= cutoff_length) || (~dorunlength && cutoff_times(i) >= cutoff_time + 0.0001)
%             %=== Even a run with a higher cutoff timed out => timeout
%             censored = 1;
%             runtime = cutoff_time;
%             runlength = cutoff_length;
%             solved = 0;
%             best_sol = best_sols(i);
%             return
%         end
    else
        %=== Successful run.
        if (dorunlength && runlengths(i) <= cutoff_length) || (~dorunlength && runtimes(i) <= cutoff_time + 0.0001)
            %=== Shorter run already successful, return its runtime.
            runtime = runtimes(i);
            runlength = runlengths(i);
            censored = 0;
            best_sol = best_sols(i);
            return
        else
% Same as above: this case is unsafe for optimizing sol. qual., so rather redo it.            
%             %=== Took longer than we have -> timeout
%             censored = 1;
%             runtime = cutoff_time;
%             runlength = cutoff_length;
%             solved = 0;
%             best_sol = best_sols(i);
%             return
        end
    end
end
runtime = [];
runlength = [];
censored = [];
best_sol = [];
solved = [];