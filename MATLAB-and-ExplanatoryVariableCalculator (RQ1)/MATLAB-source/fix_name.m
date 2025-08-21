function fixed = fix_name(name)
pairs = {};
pairs{end+1} = {'IBM-ALL-minisat', '\\minisat{}-\\IBM{}'};
pairs{end+1} = {'SWV-minisat', '\\minisat{}-\\SWV{}'};
pairs{end+1} = {'IBM-SWV-minisat', '\\minisat-\\SWVIBM{}'};
pairs{end+1} = {'INDU-minisat', '\\minisat{}-\\INDU{}'};
pairs{end+1} = {'HAND-minisat', '\\minisat{}-\\HAND{}'};
pairs{end+1} = {'RAND-minisat', '\\minisat{}-\\RAND{}'};
pairs{end+1} = {'INDU-HAND-RAND-minisat', '\\minisat{}-\\INDUHANDRAND{}'};
pairs{end+1} = {'IBM-ALL-cryptominisat', '\\cryptominisat-\\IBM{}'};
pairs{end+1} = {'SWV-cryptominisat', '\\cryptominisat-\\SWV{}'};
pairs{end+1} = {'IBM-SWV-cryptominisat', '\\cryptominisat-\\SWVIBM{}'};
pairs{end+1} = {'INDU-cryptominisat', '\\cryptominisat-\\INDU{}'};
pairs{end+1} = {'IBM-ALL-spear', '\\spear-\\IBM{}'};
pairs{end+1} = {'SWV-spear', '\\spear-\\SWV{}'};
pairs{end+1} = {'IBM-SWV-spear', '\\spear-\\SWVIBM{}'};
pairs{end+1} = {'INDU-spear', '\\spear-\\INDU{}'};
pairs{end+1} = {'RAND-SAT-saps', '\\saps{}-\\RANDSAT{}'};
pairs{end+1} = {'RAND-SAT-tnm', '\\tnm{}-\\RANDSAT{}'};
pairs{end+1} = {'BIGMIX-cplex', '\\cplex{}-\\BIGMIX{}'};
pairs{end+1} = {'BIGMIX-gurobi', '\\gurobi{}-\\BIGMIX{}'};
pairs{end+1} = {'BIGMIX-scip', '\\scip{}-\\BIGMIX{}'};
pairs{end+1} = {'BIGMIX-lpsolve', '\\lpsolve{}-\\BIGMIX{}'};
pairs{end+1} = {'CORLAT-cplex', '\\cplex{}-\\CORLAT{}'};
pairs{end+1} = {'CORLAT-gurobi', '\\gurobi{}-\\CORLAT{}'};
pairs{end+1} = {'CORLAT-scip', '\\scip{}-\\CORLAT{}'};
pairs{end+1} = {'CORLAT-lpsolve', '\\lpsolve{}-\\CORLAT{}'};
pairs{end+1} = {'RCW-cplex', '\\cplex{}-\\RCW{}'};
pairs{end+1} = {'REG-cplex', '\\cplex{}-\\REG{}'};
pairs{end+1} = {'CORLAT-REG-cplex', '\\cplex{}-\\CORLATREG{}'};
pairs{end+1} = {'CORLAT-REG-RCW-cplex', '\\cplex{}-\\CORLATREGRCW{}'};
pairs{end+1} = {'PORTGEN-lkh202', '\\lkh{}-\\PORTGEN{}'};
pairs{end+1} = {'PORTCGEN-lkh202', '\\lkh{}-\\PORTCGEN{}'};
pairs{end+1} = {'PORTGEN-PORTCGEN-lkh202', '\\lkh{}-\\PORTGENPORTCGEN{}'};
pairs{end+1} = {'TSPLIB_LKH_default_both_ok', '\\lkh{}-\\TSPLIB{}'};
pairs{end+1} = {'PORTGEN-concorde', '\\concorde{}-\\PORTGEN{}'};
pairs{end+1} = {'PORTCGEN-concorde', '\\concorde{}-\\PORTCGEN{}'};
pairs{end+1} = {'PORTGEN-PORTCGEN-concorde', '\\concorde{}-\\PORTGENPORTCGEN{}'};
pairs{end+1} = {'TSPLIB_concorde_default_both_ok', '\\concorde{}-\\TSPLIB{}'};


pairs{end+1} = {'SPEAR-ibm-swv-al', '\\spear{}-\\SWVIBM{}'};
pairs{end+1} = {'SPEAR-swv-al', '\\spear{}-\\SWV{}'};
pairs{end+1} = {'SPEAR-ibm-al', '\\spear{}-\\IBM{}'};
pairs{end+1} = {'CPLEX12-cat-BIGMIX', '\\cplex{}-\\BIGMIX{}'};
pairs{end+1} = {'CPLEX12-cat-CORLAT', '\\cplex{}-\\CORLAT{}'};
pairs{end+1} = {'CPLEX12-cat-REG', '\\cplex{}-\\REG{}'};
pairs{end+1} = {'CPLEX12-cat-RCW', '\\cplex{}-\\RCW{}'};
pairs{end+1} = {'CPLEX12-cat-CORLAT-REG', '\\cplex{}-\\CORLATREG{}'};
pairs{end+1} = {'CPLEX12-cat-CORLAT-REG-RCW', '\\cplex{}-\\CORLATREGRCW{}'};
pairs{end+1} = {'LKH-TSPLIB', '\\lkh{}-\\TSPLIB{}'};



for i=1:length(pairs)
    if strcmp(pairs{i}{1}, name)
        name = pairs{i}{2};
    end
    if strcmp(strcat('params-',pairs{i}{1}), name)
        name = pairs{i}{2};
    end
end
fixed = name;