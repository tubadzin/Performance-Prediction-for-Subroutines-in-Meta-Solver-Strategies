/*------------------------------------------------------------------------*/
/*  File: comp_feat.cpp                                                   */
/*  Feature computation for general MIP instances, including all types    */
/*  of instances CPLEX can read.                                          */
/*------------------------------------------------------------------------*/
/*  This code started as a generalization of feature computation code     */
/*  for combinatorial auction instances written by Eugene Nudelman.       */
/*  Since then, we added lots of new features.                            */
/*  The CPLEX interfacing is based on file examples/src/lpex2.c in        */
/*  CPLEX, Version 9.1                                                    */
/*  Frank Hutter, Aug 19, 2008                                            */
/*------------------------------------------------------------------------*/
/*  Differences in features compared to Nudelman's code are due to        */
/*  normalizations that were done in that code (but not reported in       */
/*  papers), and to at least one bug in that code: there is a problem in  */
/*  the bid graph degrees in that the code counts hundreds of times too   */
/*  high degrees, resulting in edge densities > 500 (should be            */
/*  upper-bounded by 1). Also, in one run valgrind found 7373 errors from */
/*  121 contexts, which suggests there may be some problem with memory.   */
/*------------------------------------------------------------------------*/
/*  I am being sloppy with freeing memory - this                          */
/*  is just fine when used as a standalone application, but if you want   */
/*  to use it in a bigger system you'll have to fix it to prevent leaks.  */
/*------------------------------------------------------------------------*/
/*  WHAT CHANGED?
/*------------------------------------------------------------------------*/
/*  RAII wrappers for CPLEX resources**									  */
/*  Added `CplexEnv`  – owns a `CPXENVptr`, opens it in the ctor,		  */
/*  auto-closes it in the dtor via `CPXcloseCPLEX`.						  */
/*  Added `CplexProb` – owns a `CPXLPptr`, creates the problem in		  */
/*  the ctor, auto-frees it in the dtor via `CPXfreeprob`.				  */



/* Bring in the declarations for the string and character functions 
and malloc, timer, and set. */

#include "comp_feat.h"
#include <cstdio>
#include <boost/lexical_cast.hpp>
#include <ilcplex/cplex.h>
#include "cplex_raii.hpp"
using namespace std;


// min and max macros
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

/* Useful macro to swap 2 elements, using
* the standard trick do { } while(0) to
* make this one statement (otherwise it
* could be wrong to use like 'else SWAP(A,B);').
*/
#define SWAP(A,B) \
do { int tmp = A; A = B; B = tmp; } while(0)

void writefeat(string name, double value){
  featname.push_back(name);
  featvalue.push_back(value);
}

void writetimefeat(string name, double value){
  timefeatname.push_back(name);
  timefeatvalue.push_back(value);
}

void printfeat(){
  vector<string>::iterator i;
  string tmpname;
  vector<double>::iterator j;
  double tmpvalue;
  for (i = featname.begin(); i != featname.end()-1; i++) {
    tmpname=*i;
    cout << tmpname <<", ";
  } 
  tmpname=featname.back();
  cout << tmpname<< endl;
  for (j = featvalue.begin(); j != featvalue.end()-1; j++) {
    tmpvalue=*j;
    cout << tmpvalue <<", ";
  } 
  tmpvalue=featvalue.back();
  cout << tmpvalue<< endl;
}

void printtimefeat(){
  vector<string>::iterator i;
  string tmpname;
  vector<double>::iterator j;
  double tmpvalue;
  for (i = timefeatname.begin(); i != timefeatname.end()-1; i++) {
    tmpname=*i;
    cout << tmpname <<", ";
  } 
  tmpname=timefeatname.back();
  cout << tmpname<< endl;
  for (j = timefeatvalue.begin(); j != timefeatvalue.end()-1; j++) {
    tmpvalue=*j;
    cout << tmpvalue <<", ";
  } 
  tmpvalue=timefeatvalue.back();
  cout << tmpvalue<< endl;
}

void terminate(int status){
  printfeat();
  cout << "\n";
  printtimefeat();
  cout << "\n";
  exit(status);
}

string convertInt(int number){
  return boost::lexical_cast<string>( number );
}

static void
free_and_null (char **ptr),
usage         (char *progname);


int size_of_intersection(int *a, int*b, int num_a, int num_b){
  /* Compute the size of the intersection of a[0],...,a[num_a-1], and b[0],...,b[num_b-1]
  This assumes a and b are unique and sorted. */
    int i=0, j=0, jstart=0, count=0;
    while (i<num_a && j<num_b){
      if( a[i] > b[j] ){
	j++;
      } else if( a[i] < b[j] ){
	i++;
      } else {
	count++;
	i++;
	j++;
      }
    }
    return count;
}

void printMatrixDouble(const double* X,int nRows, int nCols){
  int i,j;
  for(i = 0; i < nRows; i++) {
    printf("< ");
    for(j = 0; j < nCols; j++) {
      printf("%lf ",X[i+nRows*j]);
    }
    printf(">\n");
  }
}

void printMatrixInt(int* X,int nRows, int nCols){
  int i,j;
  for(i = 0; i < nRows; i++) {
    printf("< ");
    for(j = 0; j < nCols; j++) {
      printf("%d ",X[i+nRows*j]);
    }
    printf(">\n");
  }
}

int initialize(char* filename, int debug, CPXENVptr *env, CPXLPptr *lp){
  int status;
  /* Initialize the CPLEX environment */
    *env = CPXopenCPLEX (&status);
    /* If an error occurs, the status value indicates the reason for
    failure.  A call to CPXgeterrorstring will produce the text of
    the error message.  Note that CPXopenCPLEX produces no output,
    so the only way to see the cause of the error is to use
    CPXgeterrorstring.  For other CPLEX routines, the errors will
    be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */
    if ( *env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (*env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      terminate(status);
    }
    
    /* Set node log interval. */
    status = CPXsetintparam (*env, CPX_PARAM_MIPINTERVAL, 1);
    if ( status ) {
      fprintf (stderr, 
		"Failure to set CPX_PARAM_MIPINTERVAL, error %d.\n", status);
		terminate(status);
    }
    
    /* Set level pf detail for node log. */
    status = CPXsetintparam (*env, CPX_PARAM_MIPDISPLAY, 2);
    if ( status ) {
      fprintf (stderr, 
		"Failure to set CPX_PARAM_MIPDISPLAY, error %d.\n", status);
		terminate(status);
    }
    
    /* Turn on output to the screen */
    /*   status = CPXsetintparam (*env, CPX_PARAM_SCRIND, CPX_ON);
    if ( status ) {
 fprintf (stderr, 
 "Failure to turn on screen indicator, error %d.\n", status);
 terminate(status);
}*/

    /* Create the problem, using the filename as the problem name */
    *lp = CPXcreateprob (*env, &status, filename);
    /* A returned pointer of NULL may mean that not enough memory
    was available or there was some other problem.  In the case of 
    failure, an error message will have been written to the error 
    channel from inside CPLEX.  In this example, the setting of
    the parameter CPX_PARAM_SCRIND causes the error message to
    appear on stdout.  Note that most CPLEX routines return
    an error code to indicate the reason for failure.   */
    if ( *lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      terminate(status);
    }
    
    /* Now read the file, and copy the data into the created lp */
    status = CPXreadcopyprob (*env, *lp, filename, NULL);
    if ( status ) {
      fprintf (stderr, "Failed to read and copy the problem data.\n");
      terminate(status);
    }
    double time_parse = sw.Lap();
    if (debug) fprintf(stderr, " parse time = %lfs\n", time_parse);
	return status;
}

void getConstraintTypes(CPXENVptr env, CPXLPptr lp, int cur_numrows, vector< set<int> > *constraint_sets){
  char* sense = new char[cur_numrows];
  
  /* Get sense of constraints. */
    int status = CPXgetsense (env, lp, sense, 0, cur_numrows-1);
    if ( status ) {
      fprintf (stderr, "Failed to get constraint sense. Error %d\n", status);
      terminate(status);
    }
    
    (*constraint_sets) = vector< set<int> > (4); /* 'L', 'E', 'G', 'R' */
    for(int i=0; i<cur_numrows; i++){
      switch ( sense[i] ) {
	case 'L':
	  (*constraint_sets)[0].insert(i);
	  break;
	case 'E':
	  (*constraint_sets)[1].insert(i);
	  break;
	case 'G':
	  (*constraint_sets)[2].insert(i);
	  break;
	case 'R':
	  (*constraint_sets)[3].insert(i);
	  break;			
      }
    }
}


void computeLPfeatures(CPXENVptr env_relax, CPXLPptr lp_relax, int probtype, int cur_numcols){
	/* Compute LP relaxation features */
	int *indices = new int[1];
	double *lu = new double[1];
	double *up = new double[1];
	double *lu_eps = new double[1];
	double *up_eps = new double[1];
	lu[0] = 0;
	up[0] = 1;

	char *Lchar = new char[1];
	char *Uchar = new char[1];
	char *Cchar = new char[1];
	Lchar[0] = 'L';
	Uchar[0] = 'U';
	Cchar[0] = 'C';

	int i,j;
	int num_var_indices_to_make_cont;
    int *var_indices_to_make_cont;

	double sum = 0;
	double l2 = 0;
	double max_dist = -1;
	int solnstat, solnmethod, solntype;

	double   *x     = NULL;

	double *up_bound, *low_bound;
  
	/* Get upper and lower bounds. */
    up_bound = new double[cur_numcols];
    low_bound = new double[cur_numcols];
    int status = CPXgetub(env_relax, lp_relax, up_bound, 0, cur_numcols-1);
    status = CPXgetlb(env_relax, lp_relax, low_bound, 0, cur_numcols-1);
    
    /* Get variable types in xctype. */
    char *xctype = new char[cur_numcols];
    if (probtype == 0){
      /* LPs don't have ctype information, so manually fill it in. */
	    for(i=0; i<cur_numcols; i++){
	      xctype[i] = 'C';
	    }
    } else {
		status = CPXgetctype (env_relax, lp_relax, xctype, 0, cur_numcols-1);
		if ( status ) {
			fprintf (stderr, "Failed to get the type of variables, error %d.\n", status);
			terminate(status);
		}
    }


	sw2.Start();
	double lp_avg, lp_l2_avg, lp_linf, lp_objval;
	if (probtype == CPXPROB_LP || probtype == CPXPROB_QP || probtype == CPXPROB_QCP){
		/* These problem types have only continuous variables, so zero LP slack; using lp_objval=0 here */
		lp_avg = 0;
		lp_l2_avg = 0;
		lp_linf = 0;
		lp_objval = 0;
	} else {
		int *var_indices_to_make_cont = new int[cur_numcols];
		int num_var_indices_to_make_cont = 0;

		for(i=0; i<cur_numcols; i++){
			//printf("xctype[i] = %c, low_bound[i]=%lf\n", xctype[i], low_bound[i]);
			if (xctype[i] != 'C' && xctype[i] != 'S'){ // if not (continuous or semi-continuous)
				var_indices_to_make_cont[num_var_indices_to_make_cont++] = i;
				if (xctype[i] == 'B'){
					indices[0] = i;
					/* For previously binary variables, set bounds to [0,1] */
					status = CPXchgbds (env_relax, lp_relax, 1, indices, Lchar, lu);
					status = CPXchgbds (env_relax, lp_relax, 1, indices, Uchar, up);
					low_bound[i] = lu[0];
					up_bound[i] = up[0];
				} else {
					/* Set lower bound eps lower, upper bound eps higher to avoid numerical problems. */
					indices[0] = i;
					lu_eps[0] = low_bound[i]; //-CPX_INFBOUND; //
					up_eps[0] = up_bound[i]; //CPX_INFBOUND;//
					status = CPXchgbds (env_relax, lp_relax, 1, indices, Lchar, lu_eps);
					status = CPXchgbds (env_relax, lp_relax, 1, indices, Uchar, up_eps);
				}
			}
		}


		/* Relax problem: change all variables to continuous, then solve as LP (or QP/QCP) */

		/* Set LP strategy, optimize the LP and obtain solution. */
		if (probtype == CPXPROB_MILP || probtype == CPXPROB_FIXEDMILP){
			status = CPXchgprobtype (env_relax, lp_relax, CPXPROB_LP);
			status = CPXsetintparam (env_relax, CPX_PARAM_LPMETHOD, CPX_ALG_AUTOMATIC);
			if ( status ) {
				fprintf (stderr, 
						"Failed to set LP optimization method, error %d.\n", status);
				terminate(status);
			}
			status = CPXlpopt (env_relax, lp_relax);
			if ( status ) {
				fprintf (stderr, "Failed to optimize LP.\n");
				terminate(status);
			}

		} else if (probtype == CPXPROB_MIQP || probtype == CPXPROB_FIXEDMIQP){
			status = CPXchgprobtype (env_relax, lp_relax, CPXPROB_QP);
			status = CPXsetintparam (env_relax, CPX_PARAM_QPMETHOD, CPX_ALG_AUTOMATIC);
			if ( status ) {
				fprintf (stderr, 
						"Failed to set QP optimization method, error %d.\n", status);
				terminate(status);
			}
			status = CPXqpopt (env_relax, lp_relax);
			if ( status ) {
				fprintf (stderr, "Failed to optimize QP.\n");
				terminate(status);
			}

		} else {
			status = CPXchgprobtype (env_relax, lp_relax, CPXPROB_QCP);
			status = CPXsetintparam (env_relax, CPX_PARAM_LPMETHOD, CPX_ALG_AUTOMATIC);
			if ( status ) {
				fprintf (stderr, 
						"Failed to set QCP optimization method, error %d.\n", status);
				terminate(status);
			}
			status = CPXbaropt (env_relax, lp_relax);
			if ( status ) {
				fprintf (stderr, "Failed to optimize QCP.\n");
				terminate(status);
			}

		}

	    solnstat = CPXgetstat (env_relax, lp_relax);
		if      ( solnstat == CPX_STAT_UNBOUNDED ) {
			fprintf (stderr, "Relaxed model is unbounded\n");
			terminate(status);
		}
		else if ( solnstat == CPX_STAT_INFEASIBLE ) {
			fprintf (stderr, "Relaxed model could not be shown to be feasible\n");
			terminate(status);
		}
		else if ( solnstat == CPX_STAT_INForUNBD ) {
			printf ("Relaxed model is infeasible or unbounded\n");
			terminate(status);
		} 
		else {
			status = CPXsolninfo (env_relax, lp_relax, &solnmethod, &solntype, NULL, NULL);
			if ( status ) {
				fprintf (stderr, "Failed to obtain LP solution info.\n");
				terminate(status);
			}
			//printf ("Solution status %d, solution method %d\n", solnstat, solnmethod);
		}

		if ( solntype == CPX_NO_SOLN ) {
			fprintf (stderr, "Relaxed solution not available.\n");
			//goto TERMINATE;
		}
		    
		status = CPXgetobjval (env_relax, lp_relax, &lp_objval);
		if ( status ) {
			fprintf (stderr, "Failed to obtain relaxed objective value.\n"); 
			// FH: if this error occurs for an instance, one could code "failed objective" with a special value
			terminate(status);
		}
		//printf ("Objective value %lf\n", lp_objval);

		/* Retrieve solution vector */
		x = (double *) malloc (cur_numcols*sizeof(double));
		if ( x == NULL ) {
			fprintf (stderr, "No memory for solution.\n");
			terminate(status);
		}

		status = CPXgetx (env_relax, lp_relax, x, 0, cur_numcols-1);
		if ( status ) {
			fprintf (stderr, "Failed to obtain primal solution.\n");
			terminate(status);
		}

		for(j=0; j<num_var_indices_to_make_cont; j++) {
			i = var_indices_to_make_cont[j];
			double dist = min(fabs(x[i]), fabs(x[i] - 1));
			sum += dist;
			l2 += dist*dist;
			max_dist = max(max_dist, dist);
		}

		lp_avg = sum / double(num_var_indices_to_make_cont);
		lp_l2_avg = sqrt(l2/double(num_var_indices_to_make_cont));
		lp_linf = max_dist;
	}

	/* Output: LP features */
	writefeat("lp_avg", lp_avg);
	writefeat("lp_l2_avg", lp_l2_avg);
	writefeat("lp_linf", lp_linf);
	writefeat("lp_objval", lp_objval);

	double time_relax = sw2.Lap();
	writetimefeat("time_relax", time_relax);
}

void getVariableTypes(CPXENVptr env, CPXLPptr lp, int probtype, int cur_numcols, vector< set<int> > *variable_sets){
  int i;
	int num_unbounded_discrete, num_non_cont;
  double *up_bound, *low_bound;
  
  /* Get upper and lower bounds. */
    up_bound = new double[cur_numcols];
    low_bound = new double[cur_numcols];
    int status = CPXgetub(env, lp, up_bound, 0, cur_numcols-1);
    status = CPXgetlb(env, lp, low_bound, 0, cur_numcols-1);
    
    /* Get variable types in xctype. */
    char *xctype = new char[cur_numcols];
    if (probtype == 0){
      /* LPs don't have ctype information, so manually fill it in. */
	    for(i=0; i<cur_numcols; i++){
	      xctype[i] = 'C';
	    }
    } else {
      status = CPXgetctype (env, lp, xctype, 0, cur_numcols-1);
      if ( status ) {
	fprintf (stderr, "Failed to get the type of variables, error %d.\n", status);
	terminate(status);
      }
    }
    
    /* Collect integer and continuous variables into sets */
    *variable_sets = vector< set<int> > (3); /* continuous, integer&binary, all */
    int num_binary = 0;
    int num_integer = 0;
    int num_continuous = 0;
    int num_semicont = 0;
    int num_semiint = 0;
    
    for(i=0; i<cur_numcols; i++){
      switch ( xctype[i] ) {
	case 'B':
	  num_binary++;
	  break;
	case 'I':
	  num_integer++;
	  break;
	case 'C':
	  num_continuous++;
	  break;
	case 'S':
	  num_semicont++;
	  break;
	case 'N':
	  num_semiint++;
	  break;			
      }
      if (xctype[i] == 'C'){ /* not counting semi-continuous variables here! xctype[i] == 'S' */
		    (*variable_sets)[1].insert(i);
      } else {
	(*variable_sets)[0].insert(i);
      }
      (*variable_sets)[2].insert(i);
    }
    
    /* Output: basic features */
    writefeat("num_b_variables", num_binary);
    writefeat("num_i_variables", num_integer);
    writefeat("num_c_variables", num_continuous);
    writefeat("num_s_variables", num_semicont);
    writefeat("num_n_variables", num_semiint);
    writefeat("ratio_b_variables", num_binary / (cur_numcols+0.0));
    writefeat("ratio_i_variables", num_integer / (cur_numcols+0.0));
    writefeat("ratio_c_variables", num_continuous / (cur_numcols+0.0));
    writefeat("ratio_s_variables", num_semicont / (cur_numcols+0.0));
    writefeat("ratio_n_variables", num_semiint / (cur_numcols+0.0));
    
    writefeat("num_i+_variables", (*variable_sets)[0].size());
    writefeat("ratio_i+_variables", ((*variable_sets)[0].size() + 0.0) / (cur_numcols+0.0));
    
    /* Compute support sizes for integer variables and percentage of unbounded integer variables */
    vector<double> support_size ((*variable_sets)[0].size(), 0);
    if (probtype == CPXPROB_LP || probtype == CPXPROB_QP || probtype == CPXPROB_QCP){
      /* These problem types have only continuous variables, so these features are not active (they are conditional!)*/
	    num_unbounded_discrete = 0;
    } else {
      /* Percentage continuous or semi-continuous variables, and support sizes for non-continuous  variables. */
		num_non_cont = 0;
	    num_unbounded_discrete = 0;
	    set<int>::iterator it;
	    
		bool isbounded;
		for (it=(*variable_sets)[0].begin(); it!=(*variable_sets)[0].end(); it++){
			isbounded = true;
			i = *it;
			if (xctype[i] == 'B'){
				support_size[num_non_cont++] = 1+1;
			} else {
				if (xctype[i] == 'I' || xctype[i] == 'N'){
					if (up_bound[i] == CPX_INFBOUND || low_bound[i] == -CPX_INFBOUND){
						num_unbounded_discrete++;
						isbounded = false;
					} 
				}

				if (isbounded){
					if (xctype[i] == 'I'){
						support_size[num_non_cont++] = up_bound[i]-low_bound[i]+1;
					} else {
						if (xctype[i] == 'N'){ /* semi-integer */
							support_size[num_non_cont++] = up_bound[i]-low_bound[i]+1+1;
						} else {
							if (xctype[i] == 'S'){ 	/* semi-continuous: either 0 or in that interval, so "2" is the best support size I can come up with */
								support_size[num_non_cont++] = 2;
							} else {
								fprintf (stderr, "WARNING: unknown xctype at position %d: %c. Ignoring this variable in feature computation.\n", i, xctype[i]);
							}
						}
					}
				}
			}
		}
		support_size.resize(num_non_cont); // unbounded variables don't count.
    }
    writefeat("num_unbounded_disc", num_unbounded_discrete);
    if ((*variable_sets)[0].size() == 0){
      writefeat("ratio_unbounded_disc", 0);
    } else {
      writefeat("ratio_unbounded_disc", (num_unbounded_discrete+0.0)/((*variable_sets)[0].size()));
    }
    
	basicVect(support_size, noone, "support_size", standardFeatString);
}


void computeGraphFeatures(CPXENVptr env, CPXLPptr lp, set<int> variable_set, int var_set_num, bool compute_vg, bool compute_clust_coef, bool debug){
  int i,j,k_all,k;
  set<int>::iterator it;
  string var_set_numstr = convertInt(var_set_num);
  
  vector<double> A_ij_normalized(nzcnt,0.0);
  vector<double> a_normalized_varcoefs(cur_numrows,0.0);
  vector<double> vcg_constraint_degree(cur_numrows,0.0);
  vector<double> vcg_constraint_sum(cur_numrows,0.0);
  
  vector<double> obj_fun_coef_this (variable_set.size(), 0.0); 
  /* The following are really integers, but to avoid having to rewrite computation of statistics for ints, we make them doubles. */
    vector<double> vcg_var_degree(variable_set.size(), 0.0);  
    vector<double> vcg_var_sum(variable_set.size(), 0.0);
    vector<double> normed_abs_row_vals(variable_set.size(), 0.0);
    
    int *map_all_vars_to_active_vars = new int[cur_numcols];

    sw.Start();
	fprintf(stderr, "Computing VCG features - O(#nonzeros) ...\n");
    
    /* Structures for variable graph. */
    int *vg_degree = new int[variable_set.size()];
    for (k=0; k<variable_set.size(); k++) vg_degree[k] = 0;
    
    /* Compute map from all variables to active ones. VERY IMPORTANT! */
    for (k_all=0; k_all<cur_numcols; k_all++){
      map_all_vars_to_active_vars[k_all] = -1; // -1 denoting not active
    }
    int num_active=0;
    for (it=variable_set.begin(); it!=variable_set.end(); it++){
      k_all = *it;
      map_all_vars_to_active_vars[k_all] = num_active++;
    }
    
    /****************************************************************************
    * Compute variable constraint graph features.
    ****************************************************************************/
    /* Initialize counters. */
//	for (i=0; i<cur_numrows; i++) vcg_constraint_degree[i] = 0;
//	for (i=0; i<cur_numrows; i++) vcg_constraint_sum[i] = 0;
//	for (k=0; k<variable_set.size(); k++) vcg_var_degree[k] = 0;
//	for (k=0; k<variable_set.size(); k++) vcg_var_sum[k] = 0;
    
    /* Enumerate nonzero entries, get variable and constraint degrees, and compute features of A matrix */
	int count = 0;
    for(i=0; i<cur_numrows; i++){
      int for_up_bound;
      if (i==cur_numrows-1) {
	for_up_bound = nzcnt; 	/* SHOULD THIS BE -1? not sure! */ /* last row, special case for the upper limit of the loop */
      } else {
	for_up_bound = rmatbeg[i+1];
      }
      double sum_abs = 0;
      int num_nz_this_row = 0; // only consider columns in the current variable set.
      for (j=rmatbeg[i]; j<for_up_bound; j++){
	k_all = rmatind[j];
	/* If variable k_all is not an element of the set we're currently considering then skip it. */
		    k = map_all_vars_to_active_vars[k_all];
		    if (k == -1){
		      continue;
		    }
		    
		    vcg_constraint_degree.at(i)=vcg_constraint_degree.at(i)+1.0;
		    vcg_var_degree.at(k)= vcg_var_degree.at(k)+1.0; 
		    
		    /* For computing features of the rows of the A matrix */
		    double a_ik_all = rmatval[j];
		    normed_abs_row_vals[num_nz_this_row++] = fabs(a_ik_all); // absolute of nonzero entries in this row (only considering variables in the current set).
		    sum_abs += fabs(a_ik_all);
		    vcg_var_sum.at(k) = vcg_var_sum.at(k) + a_ik_all;
		    vcg_constraint_sum.at(i) =vcg_constraint_sum.at(i)+ a_ik_all;
		    
		    if (fabs(rhs_b[i] > 1e-6)){
		      /* A_ij_normalized holds constraint coefficients of normalized MIP with b = vector of ones (ignoring those rows with b=0 */
				A_ij_normalized[count++] = a_ik_all / rhs_b[i];
		    }
		    //printf(" row=%d, j=%d, k=%d, a_ik=%lf, rhs=%lf, begin=%d, end=%d,#vars in constraint=%d\n", i, j, k, a_ik, rhs_b[i], rmatbeg[i], for_up_bound , (for_up_bound - rmatbeg[i] - 1));
      }
      /* Compute variation coefficient of normalized absolute nonzero entries in this row */
	    for (j=0; j<num_nz_this_row; j++){
	      normed_abs_row_vals[j] /= sum_abs;
	    }
	    if (num_nz_this_row == 0){
			a_normalized_varcoefs[i] = 0.0;
	    } else {
			normed_abs_row_vals.resize(num_nz_this_row);
			double mysum=0.0, mymean=defaultForMissingValue, mystd=defaultForMissingValue, myvarcoef=defaultForMissingValue, myinvvarcoef=defaultForMissingValue;
			computeFeatures_O_of_N(normed_abs_row_vals, noone, mysum, mymean, mystd, myvarcoef, myinvvarcoef);
		//	cout << normed_abs_row_vals.size() << " " << mystd << " " <<mymean << endl ;
                        if (isnan(myvarcoef)){
                            display (normed_abs_row_vals);
                              cout << mysum << " " << mystd << " " << mymean << endl;
                        }
			a_normalized_varcoefs[i] = myvarcoef;
                      }	
    }

	A_ij_normalized.resize(count);
    /* Compute the stats for all the different vectors. */
	basicVect(vcg_constraint_degree, noone, "vcg_constr_deg"+ var_set_numstr, standardFeatString); // cur_numrows
	basicVect(vcg_var_degree, noone, "vcg_var_deg"+ var_set_numstr, standardFeatString); // variable_set.size()
	basicVect(vcg_constraint_sum, noone, "vcg_constr_weight"+ var_set_numstr, "avg-varcoef"); // cur_numrows
    basicVect(vcg_var_sum, noone, "vcg_var_weight"+ var_set_numstr, "avg-varcoef"); // variable_set.size()
    basicVect(A_ij_normalized, noone, "A_ij_normalized"+ var_set_numstr, "avg-varcoef");
	basicVect(a_normalized_varcoefs, noone, "a_normalized_varcoefs"+ var_set_numstr, "avg-varcoef");
    //   display(a_normalized_varcoefs);    
    /****************************************************************************
    * Compute features using objective function coefficients. 
    ****************************************************************************/
    int num_constrained_vars = 0;
    vector<double> mean_obj_coef_per_constr(variable_set.size(), 0);
    vector<double> mean_obj_coef_per_sqrt_constr(variable_set.size(), 0);
    for (it=variable_set.begin(); it!=variable_set.end(); it++){
      k_all = *it;
      k = map_all_vars_to_active_vars[k_all];
      obj_fun_coef_this[k] = fabs(obj_fun_coef[k_all]);
      if (vcg_var_degree[k] == 0){
		continue;
      } else {
		mean_obj_coef_per_constr[num_constrained_vars] = obj_fun_coef_this[k] / vcg_var_degree[k];
		mean_obj_coef_per_sqrt_constr[num_constrained_vars++] = obj_fun_coef_this[k] / sqrt(vcg_var_degree[k]);
      }
    }
    
	basicVect(obj_fun_coef_this, noone, "obj_coefs"+ var_set_numstr, "avg-std"); // variable_set.size()
	basicVect(mean_obj_coef_per_constr, noone, "obj_coef_per_constr"+ var_set_numstr, "avg-std"); // num_constrained_vars
	basicVect(mean_obj_coef_per_sqrt_constr, noone, "obj_coef_per_sqr_constr"+ var_set_numstr, "avg-std"); // num_constrained_vars
   
    double time_VCG = sw.Lap();
    writetimefeat("time_VCG" + var_set_numstr, time_VCG);
    if (debug) fprintf(stderr, " variable clause graph took %lfs\n", time_VCG);
    
	/****************************************************************************
	* Compute variable graph features. 
	****************************************************************************/
    if (compute_vg){
		fprintf(stderr, "Computing VG features - O(#num constraints * (#num vars in constraint)^2) ...\n");
		sw2.Start();

	    vector< set<int> > vset = vector< set<int> > (variable_set.size());
		int **vg_neighbours = new int*[variable_set.size()];
	  
		/* Duplicating code to cleanly split VCG and VG features */
		for(i=0; i<cur_numrows; i++){
			int for_up_bound;
			if (i==cur_numrows-1) {
				for_up_bound = nzcnt; 	/* SHOULD THIS BE -1? not sure! */ /* last row, special case for the upper limit of the loop */
			} else {
				for_up_bound = rmatbeg[i+1];
			}
			int num_nz_this_row = 0; // only consider columns in the current variable set.
			for (j=rmatbeg[i]; j<for_up_bound; j++){
				k_all = rmatind[j];
				/* If variable k_all is not an element of the set we're currently considering then skip it. */
				k = map_all_vars_to_active_vars[k_all];
				if (k == -1){
					continue;
				}

				/* VG features */
				/* Collect neighbours of variable k_all (k w.r.t. the active set) */
				for (int j2=j+1; j2<for_up_bound; j2++){
					int k2_all = rmatind[j2];
					int k2 = map_all_vars_to_active_vars[k2_all];
					if (k2 == -1){
						continue;
					}
					vset[ k ].insert( k2 );
					vset[ k2 ].insert( k );
				}
			}
		}
	
	    /* Fill vg_degree and vg_neighbours. */
	    for(k=0; k<variable_set.size(); k++){
	      vg_degree[k] = vset[ k ].size();
	      vg_neighbours[k] = new int[vg_degree[k]];
	      int count = 0;
	      for (it=vset[ k ].begin(); it!=vset[ k ].end(); it++){
		vg_neighbours[k][count++] = *it;
	      }
	    }
	    
		vector<double> vg_deg_double (variable_set.size(), 0);
		for(i=0; i<variable_set.size(); i++){
		   vg_deg_double[i] = vg_degree[i] + 0.0;
		}
		basicVect(vg_deg_double, noone, "VG_deg" + var_set_numstr, standardFeatString); // variable_set.size()
	    // printMatrixInt(vg_degree, 1, cur_numcols);
	    
	    /* Compute edge density */
	    int n_edges = 0;
	    for(k=0; k<variable_set.size(); k++){
	      n_edges += vg_degree[k];
	    }
	    n_edges /= 2;
	    double edge_density = n_edges / (double)( (variable_set.size()*(variable_set.size()-1)) / 2.0 );
	    writefeat("edge_density" + var_set_numstr, edge_density);
	    double time_VG = sw2.Lap();
	    writetimefeat("time_VG" + var_set_numstr, time_VG);
	    if (debug) fprintf(stderr, " variable graph took %lfs\n", time_VG);
	    
	    if (compute_clust_coef){
	      /****************************************************************************
	      * Compute clustering coefficient. 
	      ****************************************************************************/
		    sw2.Start();
			vector<double> all_clust_coef (variable_set.size(), 0);
		    fprintf(stderr, "Computing clustering coefficient - O(#variables^3) ...\n");
		    /* Compute clustering coefficient and deviation */
		    for (k=0; k<variable_set.size(); k++){
		      int num_neighbour_edges = 0;
		      if( vg_degree[k] <= 1 ){
			all_clust_coef[k] = 0;
		      } else {
			for (i=0; i<vg_degree[k]; i++){
			  /* For each variable k and each of its neighbours k2': intersect N(k) and N(k2') (which are already sorted), divide by 2 */
					    int k2 = vg_neighbours[k][i];
					    num_neighbour_edges += size_of_intersection(vg_neighbours[k], vg_neighbours[k2], vg_degree[k], vg_degree[k2]);
			}
			num_neighbour_edges /= 2;
			all_clust_coef[k] = num_neighbour_edges/double(vg_degree[k]*(vg_degree[k]-1)/2);
		      }
		    }
		    //	printMatrixDouble(all_clust_coef, 1, variable_set.size());
		    
			basicVect(all_clust_coef, noone, "CC" + var_set_numstr, standardFeatString); // variable_set.size()
		    double time_CC = sw2.Lap();
		    writetimefeat("time_CC" + var_set_numstr, time_CC);
		    if (debug) fprintf(stderr, " clustering coefficient took %lfs\n", time_CC);
	    }
    }	
}


/* probing features based on log file from CPLEX */
void prob_cplex(CPXENVptr env, CPXLPptr lp){
	int status=0;
	FILE* fpout  = NULL;
  CPXCHANNELptr cpxresults, cpxwarning, cpxerror, cpxlog;
  char* probfile;
  probfile=new char[256];
  //      char* probfile2;
  //      probfile2=new char[256];
  double objval=hugeValue;
  double finalgap=hugeValue;
  double mip_obj=hugeValue;
  int itcnt ;
  sw.Start();
  
  status = CPXsetdblparam (env, CPX_PARAM_TILIM, 5);
  //		status = CPXsetintparam (env, CPX_PARAM_ITLIM, 3000);
  
  sprintf(probfile, "%s", P_tmpdir);
  strcat(probfile, "/outputXXXXXX");
  mkstemp(probfile);
  fpout = fopen (probfile, "w");

  //status = CPXgetchannels (env, &cpxresults, &cpxwarning, &cpxerror, &cpxlog);
	status = CPXsetlogfilename(env, probfile, 0);
  status = CPXmipopt (env, lp);
  
  itcnt = CPXgetmipitcnt (env, lp);
  // printf("\nitcnt =%d\n",  itcnt);
  status = CPXgetmipobjval (env, lp, &mip_obj);
  
  if ( status ) {
    if (status == 1217 ){
      printf("\nmip_obj does not exist with code %d\n", status);
      mip_obj=hugeValue;
    } else {
      fprintf (stderr, "CPXgetmipobjval failed, error code %d.\n", status);
      terminate(status);
    }
  } 
  else {
    // printf("\nmip_obj = %lf\n", mip_obj);
  }
  
  status = CPXgetbestobjval (env, lp, &objval);
  if ( status ) {
    fprintf (stderr, "CPXgetbestobjval failed, error code %d.\n", status);
    terminate(status);
  }
  
 //  cout << objval << " " << mip_obj << endl; 
  if (mip_obj<hugeValue/2){
    finalgap= fabs(objval - mip_obj) / (1e-10 + fabs(objval));
    // printf("\nfinalgap = %lf\n", finalgap);
  }
  
  if (finalgap >1.0){
    writefeat("mipgap", 1.0);
  } else {
    writefeat("mipgap", finalgap);
  }
  
  int nodecnt=CPXgetnodecnt (env, lp);
  writefeat("nodecnt", nodecnt+0.0);

	int clqcnt = 0;
	CPXgetnumcuts(env, lp, CPX_CUT_CLIQUE, &clqcnt);
	writefeat("clqcnt", clqcnt+0.0);

	int covcnt = 0;
	CPXgetnumcuts(env, lp, CPX_CUT_COVER, &covcnt);
	writefeat("covcnt", covcnt+0.0);
  writetimefeat("cplex_prob_time", sw.Lap());
  
  status = CPXflushstdchannels (env);
	fclose (fpout);
 //  printf ("Finished writing to file.\n");
  
  processfile(probfile);
  remove(probfile);
}         

void processfile(char* filename){
  // write stuff into file lin
  ifstream fin; 
  string mystring;
  string targetstr;
  
  
  fin.open(filename);
  vector<string> pretimeVec;
  vector<double> pretimeDV;
  vector<string> rowVec;
  vector<string> colVec;
  vector<string> nonzeroVec;
  vector<string> cliqueMemberVec;
  vector<string> cliqueCutVec;
  vector<string> impliedBoundCuts;
  vector<string> flowCuts;
  vector<string> gomoryFractionalCuts;
  vector<string> mixedIntegerRoundingCuts;
  
  vector<string> iterVec;
  vector<string> dualObjVec;
  vector<string> reltimeVec;
  int pos;
  int col1;
  int col2;
  int col3;
  int col4;
  int col5;
  int col6; 
  int col0;
  int col7;
  int col8;
  vector<string> insolutionVec;
  vector<double> insolutionDV;
  vector<string> nodeVec;
  vector<string> nodenewVec;
  vector<string> nodeleftVec;
  vector<double> nodeDV;
  vector<double> nodenewDV;
  vector<double> nodeleftDV;
  vector<string> objVec;
  vector<string> iinfVec;
  vector<double> iinfDV;
  vector<string> bestinVec;
  vector<double> bestinDV;        
  vector<string> bestnodeVec;
  vector<double> bestnodeDV;
  vector<string> itcntVec;
  vector<double> itcntDV;
  vector<string> gapVec;
  vector<double> gapDV;
  vector<string> cutsVec;
  vector<double> cutsDV;
  vector<string> numcutsVec;
  vector<double> numcutsDV;
  vector<string> integralVec;
  vector<double> integralDV;
  vector<string> nodeobjVec;
  vector<double> nodeobjDV;
  
/* TODO: if string not found for the cuts, set it to zero. Output! */
  
  while (! fin.eof() ){
    getline(fin, mystring);
//    cout << "The string is " << mystring << endl;
    targetstr="Presolve time =";
    getinfo(mystring, targetstr, pretimeVec);
    string targetstr="Reduced MIP has ";
    getinfo(mystring, targetstr, rowVec);
    targetstr="rows, ";
    getinfo(mystring, targetstr, colVec);
    targetstr="columns, and ";
    getinfo(mystring, targetstr, nonzeroVec);
    targetstr="Clique table members:";
    getinfo(mystring, targetstr, cliqueMemberVec);
    targetstr="Clique cuts applied:";
    getinfo(mystring, targetstr, cliqueCutVec);
	targetstr = "Implied bound cuts applied:";
	getinfo(mystring, targetstr, impliedBoundCuts);
	targetstr = "Flow cuts applied:";
	getinfo(mystring, targetstr, flowCuts);
	targetstr = "Gomory fractional cuts applied:";
	getinfo(mystring, targetstr, gomoryFractionalCuts);
	targetstr = "Mixed integer rounding cuts applied:";
	getinfo(mystring, targetstr, mixedIntegerRoundingCuts);
    targetstr="Iteration: ";
    getinfo(mystring, targetstr, iterVec);
    targetstr="Dual objective     =   ";
    getinfo(mystring, targetstr, dualObjVec);
    targetstr="Root relaxation solution time =  ";
    getinfo(mystring, targetstr, reltimeVec);
    //now, we need to deal with the table. (note: they may have two or more of tables)
    targetstr="   Node ";
    pos = mystring.find(targetstr);
    if (pos>=0){
      //      cout << "find the table \n";
      col0=0;
      targetstr="   Node";
      col1=mystring.find(targetstr)+targetstr.size();
      targetstr=" Left";
      col2=mystring.find(targetstr)+targetstr.size();   
      targetstr=" Objective";
      col3=mystring.find(targetstr)+targetstr.size();  
      targetstr=" IInf";
      col4=mystring.find(targetstr)+targetstr.size();  
      targetstr=" Best Integer";
      col5=mystring.find(targetstr)+targetstr.size();  
      targetstr=" Best Node";
      col6=mystring.find(targetstr)+targetstr.size();  
      targetstr=" ItCnt";
      col7=mystring.find(targetstr)+targetstr.size();  
      targetstr=" Gap";
      col8=mystring.find(targetstr)+targetstr.size();    
      getline(fin, mystring);
//      cout << "The string is " << mystring << endl;
      int take=1;
      while (take>0){
	getline(fin, mystring);
//	cout << "The string is " << mystring << endl;
	pos=mystring.find("                                                                ");
	int pos1=mystring.find("CPLEX Error");
	if (pos>=0 || mystring.empty() || pos1 <0){
	  take=0;
	  break;
	}
	pos=mystring.find("Elapsed");
	
	if (pos<0){ // This line begins with "Elapsed" or with "CPLEX ERROR" -> ignore it.
	  getinfosubstr(mystring, 0, 1,  insolutionVec);
	  getinfosubstr(mystring, 1, col1-1,  nodeVec);
	  getinfosubstr(mystring, col1, 1,  nodenewVec);
	  getinfosubstr(mystring, col1+1, col2-col1,  nodeleftVec);
	  getinfosubstr(mystring, col2+1, col3-col2,  objVec);
	  getinfosubstr(mystring, col3+1, col4-col3,  iinfVec);
	  getinfosubstr(mystring, col4+1, col5-col4,  bestinVec);
	  getinfosubstr(mystring, col5+1, col6-col5,  bestnodeVec);
	  getinfosubstr(mystring, col6+1, col7-col6,  itcntVec);
	  getinfosubstr(mystring, col7+1, col8-col7-1,  gapVec);
	}                       
      }
    }
    
  }
  fin.close();
  
  int type1 = 1;
  int type2 = 2;
  double bigvalue=1e100;
  double smallvalue = 0.0;
  
  str2double(itcntVec, itcntDV, type1, smallvalue, smallvalue);
  basicVect(itcntDV, noone, "itcnt", "max");
  
  str2double(insolutionVec, insolutionDV, type2, smallvalue, smallvalue);
  basicVect(insolutionDV, noone, "numnewsolution", "sum");

/* dominated by the more accurate reading of nodecnt   
  str2double(nodeVec, nodeDV, type1, smallvalue, smallvalue);
  basicVect(nodeDV, noone, "nodeprob", "max"); */
  
  str2double(nodenewVec, nodenewDV, type2, smallvalue, smallvalue);
  basicVect(nodenewDV, noone, "newin", "sum");
  
  str2double(nodeleftVec, nodeleftDV, type1, smallvalue, smallvalue);
  basicVect(nodeleftDV, smallvalue, "nodeleft", "avg-varcoef");
  
  splitoneVect(objVec, nodeobjVec, integralVec);
  
  str2double(nodeobjVec, nodeobjDV, type1, hugeValue, hugeValue);
  
  vector <double> diffObj;
  diffVect(nodeobjDV, diffObj, hugeValue);
  basicVect(diffObj, hugeValue, "diffObj", standardFeatString);

    str2double(integralVec, integralDV, type2, smallvalue, smallvalue);
    double counter;
    countitemVect(integralDV, 1.0, counter);
    writefeat("numfeas", counter);
    /* TODO: divide by steps; many times. */

	    str2double(iinfVec, iinfDV, type2, smallvalue, smallvalue);
	    basicVect(iinfDV, smallvalue, "iinf", standardFeatString);
	    
	    str2double(bestinVec, bestinDV, type1, hugeValue, hugeValue);
		vector <double> diffBestInt;
		diffVect(bestinDV, diffBestInt, hugeValue);
	    basicVect(diffBestInt, hugeValue, "diffBestInt", standardFeatString);
	    
	    splitVect(bestnodeVec, cutsVec, numcutsVec);
	    str2double(bestnodeVec, bestnodeDV, 1, hugeValue, hugeValue);
		vector <double> diffBestObjUp;
		diffVect(bestnodeDV, diffBestObjUp, hugeValue);
	    basicVect(diffBestObjUp, hugeValue, "diffBestObjUp", standardFeatString);
	    
	    str2double(numcutsVec, numcutsDV, 1, hugeValue, hugeValue);
	    basicVect(numcutsDV, hugeValue, "numcuts", "sum"); // TODO: maybe normalize by itcnt ?
	    
	    str2double(gapVec, gapDV, 1, hugeValue, hugeValue);
		vector <double> diffGap;
		diffVect(gapDV, diffGap, hugeValue);
	    basicVect(diffGap, hugeValue, "diffGap", standardFeatString);
	    
	    writefeat("pre_t", readfirst(pretimeVec, defaultForMissingValue));
	    writefeat("rel_t", readfirst(reltimeVec, defaultForMissingValue));
	    writefeat("new_row", readfirst(rowVec, defaultForMissingValue));
	    writefeat("new_col", readfirst(colVec, defaultForMissingValue));
		/* TODO: reduction factors of #constraints and #variables */
	    writefeat("new_nonzero", readfirst(nonzeroVec, defaultForMissingValue));
	    writefeat("clique_table", readfirst(cliqueMemberVec, defaultForMissingValue));       

	    writefeat("cliqueCuts", readfirst(cliqueCutVec, 0));       
	    writefeat("impliedBoundCuts", readfirst(impliedBoundCuts, 0));       
	    writefeat("flowCuts", readfirst(flowCuts, 0));       
	    writefeat("mixedIntegerRoundingCuts", readfirst(mixedIntegerRoundingCuts, 0));       
	    writefeat("gomoryFractionalCuts", readfirst(gomoryFractionalCuts, 0));       
}

// difference between two item
void diffVect(vector<double> dv, vector<double> &dv2, double notcount){
  double tmp, tmp2;
  vector<double>::iterator i;
  if (dv.size()>1){
    for (i = dv.begin(); i != dv.end()-1; i++) {
      tmp = *i;
      tmp2 = *(i+1);
      if ((fabs(tmp-notcount)<1e-10) || (fabs(tmp2-notcount)<1e-10))
	dv2.push_back(notcount);
      else
	dv2.push_back(tmp2-tmp);
    }
  }
}

void computeFeatures_O_of_N(vector<double> array, double notcount, double &mysum, double &mymean,  double &mystd, double &myvarcoef, double &myinvvarcoef){
  double tmp;  
  int index=0;
  vector<double> foo;     
  vector<double>::iterator i;
 if (abs(notcount-noone)>1e-10){
  for (i = array.begin(); i != array.end(); i++) {
    tmp=*i;
    if (abs(tmp-notcount)>1e-10){
      foo.push_back(tmp);
    }
  }
}
 else
    foo=array;

  double running_mean=mymean;
  double running_var_times_N=mystd;
  index=foo.size();
  int k = 2;
  if (index > 0){
	 mysum = foo[0];
     running_mean = foo[0];
	 running_var_times_N = 0;
	for (i = foo.begin()+1; i != foo.end(); i++) {
      tmp=*i;
	  double oldmean = running_mean;
	  running_mean += (tmp-running_mean)/k;
	  running_var_times_N += (tmp-oldmean)*(tmp-running_mean); // Sk = Sk-1 + (xk - Mk-1)*(xk - Mk).
      mysum=mysum+tmp;
	  k++;
    }
    mymean= mysum/index;
    mystd=sqrt(running_var_times_N/index);
  //  cout << mystd <<" " << (mysqrsum-mymean*index*mymean)/mysqrsum <<" "<<endl;
    //pay attention to num. problem.
    //    mystd=max(mystd, 0.0001);

	double mean_for_varcoef = 0;
	if (fabs(mymean) < 1e-10){
		mean_for_varcoef = 1e-10;
	} else {
		mean_for_varcoef = mymean;
	}
	myvarcoef = mystd/mean_for_varcoef;

	double std_for_invvarcoef = 0;
	if (fabs(mystd) < 1e-10){
		std_for_invvarcoef = 1e-10;
	} else {
		std_for_invvarcoef = mystd;
	}
	myinvvarcoef = mymean/std_for_invvarcoef;
  }
}

void computeFeatures_O_of_NlogN(vector<double> array, double notcount, double &mymax, double &mymin, double &myq10, double &myq25, double &mymedian, double &myq75, double &myq90 ){
  double tmp;  
  int index=0;
  vector<double> foo;     
  vector<double>::iterator i;
 if (abs(notcount-noone)>1e-10){
  for (i = array.begin(); i != array.end(); i++) {
    tmp=*i;
   if (abs(tmp-notcount)>1e-10){
      foo.push_back(tmp);
    }
  }
 }
 else
    foo=array;
 
  index=foo.size();
  if( index > 0){
	  sort(foo.begin(), foo.end());  
	  mymin=foo[0];
		myq10=foo[index/10];
		myq25=foo[index/4];
		mymedian=foo[index/2];
		myq75=foo[index*3/4];
		myq90=foo[index*9/10];
		mymax=foo[index-1];
	}
}



//write a set of feat 
void writeMfeature(string name, string flist, double mysum, double mymean, double mymax, double mymin, double mystd, double myq10, double myq25, double mymedian, double myq75, double myq90, double myvarcoef, double myinvvarcoef){
  string tmp, onename;
  vector<string> feat;
  vector<string>::iterator i;
  int j,pos;
  double foo;
  
  for (j=0; j<=20; i++){
    pos=flist.find("-");
    if (pos<0){
      feat.push_back(flist); //donot forget the last one DOG!mak
      break;
    }
    else {
      feat.push_back(flist.substr(0, pos));
      flist=flist.substr(pos+1);
    }
  }
  
  double fvalue;
  for (i = feat.begin(); i != feat.end(); i++) {
    tmp=*i;  
     
    if (tmp.find("avg")!=string::npos){
      onename = name + "_" + tmp;      
      fvalue= mymean;
    }
	if (tmp.find("varcoef")!=string::npos){
      onename = name + "_" + tmp;      
       fvalue=  myvarcoef;
    }
	if (tmp.find("invvarcoef")!=string::npos){
      onename = name + "_" + tmp;      
       fvalue= myinvvarcoef;
    }
    if (tmp.find("sum")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue= mysum;
    } 
    if (tmp.find("std")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue= mystd;
    } 
    if (tmp.find("min")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue=  mymin;
    } 
    if (tmp.find("max")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue=  mymax;
    } 
    if (tmp.find("median")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue= mymedian;
    } 
    if (tmp.find("q10")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue= myq10;
    } 
    if (tmp.find("q25")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue=  myq25;
    } 
    if (tmp.find("q75")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue=  myq75;
    }  
    if (tmp.find("q90")!=string::npos){
      onename = name + "_" + tmp;            
       fvalue=  myq90;
    } 
    if (tmp.find("q75dq25")!=string::npos){
      onename = name + "_" + tmp;            
      if (fabs(myq75-defaultForMissingValue)<1e-10 || fabs(myq25-defaultForMissingValue)<1e-10)
         foo=defaultForMissingValue;
      else{
        if (myq75< 1e-6) {
            foo=0;
        } else {
			if( myq25< 1e-6) {
				foo = defaultForMissingValue;
			} else {
				foo=myq75/myq25;
			}
		}
      }  
       fvalue= foo; 
   }
  if (tmp.find("q75mq25")!=string::npos){
      onename = name + "_" + tmp;            
	  if (fabs(myq75-defaultForMissingValue)<1e-10 || fabs(myq25-defaultForMissingValue)<1e-10)
         foo=defaultForMissingValue;
      else{
        foo=myq75-myq25;
      }  
       fvalue=  foo;
    }    
  if (tmp.find("q90mq10")!=string::npos){
      onename = name + "_" + tmp;            
    	if (fabs(myq90-defaultForMissingValue)<1e-10 || fabs(myq10-defaultForMissingValue)<1e-10)
         foo=defaultForMissingValue;
      else{
        foo=myq90-myq10;
      }  
       fvalue=  foo;
    }    
  if (tmp.find("maxmmin")!=string::npos){
      onename = name + "_" + tmp;            
    	if (fabs(mymax-defaultForMissingValue)<1e-10 || fabs(mymin-defaultForMissingValue)<1e-10)
         foo=defaultForMissingValue;
      else{
        foo=mymax-mymin;
      }  
      fvalue=  foo;
    }  
writefeat(onename, fvalue);  
  }
    
}

//vector basic information
void basicVect(vector<double> &dv, double notcount, string feat_name, string which_statistics){

  double mymin=defaultForMissingValue;
  double mymax=defaultForMissingValue;
  double mysum=0.0;
  double mysqrsum=0.0;
  double mymean=defaultForMissingValue;
  double mystd=defaultForMissingValue;
  double myq10=defaultForMissingValue;
  double myq25=defaultForMissingValue;
  double myq75=defaultForMissingValue;
  double myq90=defaultForMissingValue;
  double mymedian=defaultForMissingValue;
  double myvarcoef=defaultForMissingValue;
  double myinvvarcoef=defaultForMissingValue;
  if  ((which_statistics.find("q")!=string::npos) || (which_statistics.find("min")!=string::npos) || (which_statistics.find("max")!=string::npos) ||(which_statistics.find("median")!=string::npos) )
      computeFeatures_O_of_NlogN(dv, notcount, mymax, mymin, myq10, myq25, mymedian, myq75, myq90);

  computeFeatures_O_of_N(dv, notcount, mysum, mymean, mystd, myvarcoef, myinvvarcoef);

  writeMfeature(feat_name, which_statistics, mysum, mymean, mymax, mymin, mystd, myq10, myq25, mymedian, myq75, myq90, myvarcoef, myinvvarcoef);
  
}

void countitemVect(vector<double> dv, double item, double & numitem){
  double tmp;
  numitem=0.0;
  vector<double>::iterator i;
  for (i = dv.begin(); i != dv.end(); i++) {
    tmp = *i;
    if (fabs(tmp-item)<1e-10){
      numitem++;
    }
  }
  
}



//split vector string for cuts/best node
void splitVect(vector<string> &sv, vector<string> &sv1, vector<string> &sv2){	
  vector<string>::iterator i;
  string s1, s2, tmp;
  int pos;
  
  for (i = sv.begin(); i != sv.end(); i++) {
    tmp=*i;
    if (tmp.find(":")==string::npos){
      sv1.push_back("");
      sv2.push_back("");
      
    }
    else {
      pos=tmp.find(":");
      sv1.push_back(tmp.substr(0, pos+1));
      sv2.push_back(tmp.substr(pos+1));
      *i="";
      //                   cout << sv1.back() << " " << sv2.back() << endl;
    }
  }
}

//split vector string for objection since it may contain integral etc.
void splitoneVect(vector<string> &sv, vector<string> &sv1, vector<string> &sv2){	
  vector<string>::iterator i;
  string tmp;
  
  for (i = sv.begin(); i != sv.end(); i++) {
    tmp=*i;
    if (tmp.find(".")==string::npos){
      sv1.push_back("");
      sv2.push_back(tmp);
    }
    else {
      sv1.push_back(tmp);
      sv2.push_back("");
    }
  }
}

// transfer string to doulble with different refill option
// 1: refill with previous value
// 2: refill with bigvalue

void str2double(vector<string> &sv, vector<double> &dv, int refill, double firstline, double bigvalue){
  vector<string>::iterator i;
  double tmp;
  double foo;
  string stmp;
  if (refill == 1){    
    for (i = sv.begin(); i != sv.end(); i++) {
      string stmp=*i;
      strip(stmp);
      if (stmp.empty()){	
	if (i !=sv.begin()){
	  foo=dv.back();
	  dv.push_back(foo);}
	  else
	    dv.push_back(firstline);
      }
      else {
	tmp =  str_to_f(*i);
	dv.push_back(tmp);
      }
    }
  }
  if (refill ==2){
    for (i = sv.begin(); i != sv.end(); i++) {
      stmp=*i;
      strip(stmp);
      
      if (stmp.empty()){
	if (i !=sv.begin()){
	  dv.push_back(bigvalue);}
	  else
	    dv.push_back(firstline);
      }
      else {
	tmp =  str_to_f(*i);
	dv.push_back(tmp);
      }
    }
  }
  if (refill >2)
    cout << "Error: Refill must be 1 or 2! \n";
  
}

void getinfosubstr(string mystring, int ss, int ee,  vector<string> &sv){
  string tmpstr;
  string newstr;
  newstr=mystring.substr(ss, ee);
  sv.push_back(newstr);
}

void getinfo(string mystring, string targetstr, vector<string> &sv){
  int pos1;
  pos1 = mystring.find (targetstr);
  string newstr;
  string tmpstr;
  if (pos1>=0) {
    newstr=mystring.substr(pos1+targetstr.size());
    stringstream ss(newstr);
    ss >> tmpstr;
    // strip(tmpstr);
    sv.push_back(tmpstr);
  }
}

void display(vector<double> &sv){
  vector<double>::iterator i;
  for (i = sv.begin(); i != sv.end(); i++) {
    cout << *i << endl;
  }
  cout << "=================================" << endl;
}

void strip(string &s) {
  string::iterator i = s.begin();
  while (isspace(*i)) {
    s.erase(i);
    i = s.begin();
  }
  string::reverse_iterator j = s.rbegin();
  while (isspace(*i)) {
    s.resize(s.length()-1);
    j = s.rbegin();
  }
}


double readfirst(vector<string> dv, double valueForEmpty){
  string tmp;
  double myvalue=valueForEmpty;
  vector<string>::iterator i;
  for (i = dv.begin(); i != dv.end(); i++) {
    tmp = *i;
    myvalue=str_to_f(tmp);
    return myvalue;
  }
  return myvalue;
}

// transfer a string to double
double str_to_f(string newstr){
  double myvalue ; 
  string tmp;
  istringstream in(newstr);
  double one=1.0;
  double two=2.0;
  double three=3.0;
  double four=4.0;
  double unknown=9.0;
  // in >> tmp;
  
  if (newstr.find("*")!=string::npos){
    return one;
  }
  if (newstr.find("+")!=string::npos)
      return one;
  
  
  if (newstr.find("Cuts")!=string::npos)
      return one;
  //if (newstr.find("Impl")!=string::npos)
      //    return two;
      
      if (newstr.find("Cliques")!=string::npos)
      return two;
      
      if (newstr.find("integral")!=string::npos)
      return one;
      if (newstr.find(":")!=string::npos)
	return unknown;
      
      
      in >> myvalue;
      return myvalue;
}


// Now it is time for main

int main (int argc, char *argv[]){
  /* Declare variables. */
	CplexEnv  env;              // opens & auto-closes
	CplexEnv  env_relax;        // separate env for LP relaxation
	CplexProb lp(env,        argv[1]);   // reads file inside ctor later
	CplexProb lp_relax(env_relax, argv[1]);

    int probtype;
    
    vector< set<int> > variable_sets;
    vector< set<int> > constraint_sets;
    
    int           status = 0;
    int           i,j,k;
    int nzcnt_returned, surplus;
    
    bool compute_vg = false;
    bool compute_clust_coef = true;
    bool compute_relax = true;
    bool debug = false;
    
    sw.Start();
    
    /* Check the command line arguments */
    if ( argc != 2 ) {
      usage (argv[0]);
      terminate(status);
    }
    
    /* Initialize the problem. */

    /* Get problem type. */
    probtype = CPXgetprobtype (env, lp);
    writefeat("probtype", probtype);
    
    /* FH TODO: we're outputting basic information here. This information
    should instead be computed and outputted for all variable subsets. */

    /* Ask CPLEX for the problem size.  cur_numrows and cur_numcols store the 
    current number of rows and columns, respectively.  */
    cur_numcols = CPXgetnumcols (env, lp);
    cur_numrows = CPXgetnumrows (env, lp);
    nzcnt = CPXgetnumnz(env, lp);   /* this does *not* include the nonzero entries for quadratic constraints */
    writefeat("n_vars", cur_numcols);
    writefeat("n_constr", cur_numrows);
    writefeat("n_nzcnt", nzcnt);
    
    fprintf(stderr, "#vars=%d, #constraints=%d, #nonzeros=%d\n", cur_numcols, cur_numrows, nzcnt);
    
    /* Compute features for quadratically-constrained problems */
    int num_quad = CPXgetnumquad (env, lp); // number of variables that have quadratic objective coefficients
    int num_q_constr = CPXgetnumqconstrs (env, lp);
    int num_qpnz = CPXgetnumqpnz (env, lp);
    writefeat("nq_vars", num_quad);
    writefeat("nq_constr", num_q_constr);
    writefeat("nq_nzcnt", num_qpnz);
    
    /****************************************************************************
    * Get information about A and b.
    ****************************************************************************/
    /* Get objective function coefficients. */
    obj_fun_coef = new double[cur_numcols];
    status = CPXgetobj (env, lp, obj_fun_coef, 0, cur_numcols-1);
    if ( status ) {
      fprintf (stderr, "Failed to get objective function coefficients. Error %d\n", status);
      terminate(status);
    }
    
    /* Get right hand side, the b array. */
    rhs_b = new double[cur_numrows];
    status = CPXgetrhs (env, lp, rhs_b, 0, cur_numrows-1);
    if ( status ) {
      fprintf (stderr, "Failed to get right hand size, the b array. Error %d\n", status);
      terminate(status);
    }
    
    /* Get rmatbed, rmatind, and rmatval. */
    rmatbeg = new int[cur_numrows];
    rmatind = new int[nzcnt];
    rmatval = new double[nzcnt];
    status = CPXgetrows(env, lp, &nzcnt_returned, rmatbeg, rmatind, rmatval, nzcnt, &surplus, 0, cur_numrows-1);
    if( nzcnt_returned != nzcnt ){
      printf("nzcnt = %d, nzcnt_returned = %d\n", nzcnt, nzcnt_returned);
      fprintf(stderr, "Error in computation of nonzero entries of the constraint matrix -- value returned by CPXgetrows different than the one from CPXgetnumnz\n" );
      terminate(status);
    }
    if( surplus < 0 ){
      fprintf(stderr, "negative surplus %d, array too small\n", surplus);
      terminate(status);
    }	
    
	/* Compute LP features. */
	computeLPfeatures(env_relax, lp_relax, probtype, cur_numcols);

    /* Get variable types. */
    getVariableTypes(env, lp, probtype, cur_numcols, &variable_sets);
    
    /* Get constraint sets (<=, =, >=, 'R' */
    getConstraintTypes(env, lp, cur_numrows, &constraint_sets);
    
    
    /****************************************************************************
    * Compute basic information for each constraint type.
    ****************************************************************************/
    for (int const_set_num = 0; const_set_num<3; const_set_num++){
      string const_set_numstr = convertInt(const_set_num);
      
      /* Get RHS vector for this constraint type. */
	    vector<double> rhs_b_this(constraint_sets[const_set_num].size(), 0);
	int num_rhs_b_this = 0;
	    for (set<int>::iterator it=constraint_sets[const_set_num].begin(); it!=constraint_sets[const_set_num].end(); it++){
	      i = *it;
			rhs_b_this[num_rhs_b_this++] = rhs_b[i];
	    }
		
		basicVect(rhs_b_this, noone, "rhs_c_" + const_set_numstr, "avg-varcoef"); 
	}
    
    
    /****************************************************************************
    * Compute information for each variable type.
    ****************************************************************************/
    computeGraphFeatures(env, lp, variable_sets[0], 0, false, false, debug);
    computeGraphFeatures(env, lp, variable_sets[1], 1, false, false, debug);
    computeGraphFeatures(env, lp, variable_sets[2], 2, false, false, debug);
    
    /****************************************************************************
    * Probing features from running CPLEX for a short period of time.
    ****************************************************************************/
	prob_cplex(env, lp);       

	terminate(0);
}

/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL */
static void
free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
} /* END free_and_null */

static void
usage (char *progname)
{
  fprintf (stderr,"Usage: %s filename\n", progname);
  fprintf (stderr,"   where filename is a MIP instance file in .MPS, .LP, or .SAV format.\n");
  fprintf (stderr," Exiting...\n");
} /* END usage */
