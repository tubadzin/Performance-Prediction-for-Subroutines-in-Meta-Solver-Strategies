Note: you need CPLEX to run the feature computation code.

For detailed information, please contact 
Lin Xu, xulin730@cs.ubc.ca
Frank Hutter, hutter@cs.ubc.ca

===============================================
Short description of features:

num_b_variables ... num_n_variables: Number of variables of type (boolean, integer, continuous, semi-continuous, semi-integer)
ratio_b_variables ... ratio_n_variables: Ratio of each type of variables (summing to 1).
Define: V=all variables; C=continuous variables; O = V\C; X: current one (either V, O, or C)
num_i+_variables, ratio_i+_variables: |O| and |O|/|V|
num_unbounded_disc, ratio_unbounded_disc: x=|{v in O| unbounded on either side|, x/|O|
edge_density: #edges in VG graph/#edges in a complete graph with that many nodes

probtype: categorical feature, CPLEX problem type
n_vars, n_constr, n_nzcnt: #variables, #constraints, #nonzeros in A
nq_vars, nq_constr, nq_nzcnt: #variables with quadratic constraints; #quadratic constraints; #nonzeros in Q

Probing features (5s of CPLEX)
mipgap: MIP gap
nodecnt: # nodes visited
clqcnt, covcnt, : #clique and cover cuts applied
numfeas: #feasible solutions found
pre_t: presolve time
rel_t: relaxation time
new_row: #constraints after presolving
new_col: #variables after presolving
new_nonzero: #nonzero entries in A after presolving
clique_table: size of clique table (whatever that is ???)

Vectors: 
support_size: only for bounded: domain size for binary or integer variables, 2 for semi-cont, 1+domain size for semi-int
-> standard

vcg_constraint_degree: VCG constraint degree: (#{A_{i,j} | A_{i,j} \neq 0, x_i \in X})
-> standard

vcg_var_degree: VCG variable degree: (#{A_{i,j} | A_{i,j} \neq 0, c_j \in C})
-> standard

M = #constraints
vcg_var_sum: VCG variable sum: (sum_{j=1}^M A_{i,j})_{i \in 1, ..., |X|}
-> "avg-varcoef"

vcg_constraint_sum: VCG constraint sum: (sum_{i=1}^{|X|} A_{i,j})_{j \in 1,...,M}
-> "avg-varcoef"

A_ij_normalized: 
-> "avg-varcoef"

a_normalized_varcoefs: varcoef( (abs(A_{ij}) / sum_{i=1}^{|X| abs(A_{ij})) )_{j \in 1,...,M}
-> "avg-varcoef"

obj_coefs, mean_obj_coef_per_constr, mean_obj_coef_per_sqrt_constr: (|c_i|)_{i \in X}, (|c_i|/n_i)_{i \in X}, (|c_i|/sqrt(n_i))_{i \in X}, where n_i is # of nonzero entries in column i of A,
-> "avg-varcoef"

vg_deg_double: VG degree
-> standard

all_clust_coef: CC 
-> standard

rhs_b_this
-> {b_j}_{j in constraints with sense S}. 3 different vectors, one each for S in (<=, =, >=)
-> "avg-varcoef"



Probing features:
itcnt-max: number of iterations
numnewsolution-sum: #times that CPLEX found a new solution or incumbent
nodeprob-max: 
newin-sum: #times that CPLEX found a new incumbent by primal heuristics
nodeleft: nodes left during the trajectory, 
-> "avg-varcoef"

diffObj: improvement of objective function per output line in the log
-> standard

iinf: number of integer-infeasible variables at current node
-> standard

diffBestInt: improvement of best integer solution per output line in the log
-> standard

diffBestObjUp: improvement of upper bound per output line in the log
-> standard

numcuts-sum: total number of cuts applied

diffGap: improvement of gap per output line in the log
-> standard


NOTE: -512 encodes features that could not be computed

