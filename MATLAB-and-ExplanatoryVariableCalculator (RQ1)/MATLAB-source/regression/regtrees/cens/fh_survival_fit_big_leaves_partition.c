/* fh_random_regtreefit_logrank_big_leaves_twofeaturetypes
Builds a regression tree from the input, picking features partly at random.
*/
 
#include "mex.h"
#include <stdlib.h> 
#include <math.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* Useful macro to swap 2 elements, using
 * the standard trick do { } while(0) to
 * make this one statement (otherwise it
 * could be wrong to use like 'else SWAP(A,B);').
 */
#define SWAP(A,B) \
do { int tmp = A; A = B; B = tmp; } while(0)

#define D_SWAP(A,B) \
do { double tmp = A; A = B; B = tmp; } while(0)


/*=== Global variable for using built-in qsort since qsort_r can't be linked here.*/
double* ptr_to_double_array_for_qsort;

void quick(int *a, int min, int max) {
  if (max - min > 1) {

    int i = min;
    int j = max;
    /* indices are positive */
    int pivot = a[(i+j) >> 1];
    do {
      while(a[i] < pivot) i++;
      while(a[j] > pivot) j--;
      if (i > j) break;
      SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && a[j] == pivot) j--;
    if (min < j) quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && a[i] == pivot) i++;
    if (i < max) quick(a, i, max);

  } else if (a[min] > a[max])
    SWAP(a[min], a[max]);
}


void d_quick(double *a, int min, int max) {
  if (max - min > 1) {
    int i = min;
    int j = max;
    /* indices are positive */
    double pivot = a[(i+j) >> 1];
    do {
      while(a[i] < pivot) i++;
      while(a[j] > pivot) j--;
      if (i > j) break;
      D_SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && a[j] == pivot) j--;
    if (min < j) d_quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && a[i] == pivot) i++;
    if (i < max) d_quick(a, i, max);
  } else if (a[min] > a[max])
    D_SWAP(a[min], a[max]);
}


void dp_quick(int *a, int min, int max) {
  if (max - min > 1) {

    int i = min;
    int j = max;
    /* indices are positive */
    double pivot = ptr_to_double_array_for_qsort[a[(i+j) >> 1]];
    do {
      while(ptr_to_double_array_for_qsort[a[i]] < pivot) i++;
      while(ptr_to_double_array_for_qsort[a[j]] > pivot) j--;
      if (i > j) break;
      SWAP(a[i], a[j]);
    } while(++i <= --j);

    /* Try to reduce bad behaviours. */
    while (min < j && ptr_to_double_array_for_qsort[a[j]] == pivot) j--;
    if (min < j) dp_quick(a, min, j);

    /* Try to reduce bad behaviours. */
    while (i < max && ptr_to_double_array_for_qsort[a[i]] == pivot) i++;
    if (i < max) dp_quick(a, i, max);

  } else if (ptr_to_double_array_for_qsort[a[min]] > ptr_to_double_array_for_qsort[a[max]])
    SWAP(a[min], a[max]);
}

/* Quick sort for sorting an int array a.
   Example: quick_sort(a, size);
// Equivalent: qsort(a, size, sizeof(int), compare_ints);
*/
void quick_sort(int *a, int size) {
  if (size > 1) {
    quick(a, 0, size-1);
  }
}

/* Quick sort for sorting a double array b.
   Example: d_quick_sort(b, size);
// Equivalent: qsort(b, size, sizeof(double), compare_doubles);
*/
void d_quick_sort(double *a, int size) {
  if (size > 1) {
    d_quick(a, 0, size-1);
  }
}

/* Quick sort for sorting the indices a of a double array b.
   Example: (int array a, double arrays b, c)
	for(i=0; i<size; i++) a[i]=i;
	ptr_to_double_array_for_qsort = b;
    dp_quick_sort(a, size);
//	Equivalent: qsort(a, size, sizeof(int), compare_idxs_into_double_array);

    If we also want the double array sorted, get it in array c:
    for(i=0; i<size; i++) c[i]=b[a[i]];*/
void dp_quick_sort(int *a, int size) {
  if (size > 1) {
    dp_quick(a, 0, size-1);
  }
}

/* For randperm in C */
/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(int *array, size_t n)
{
    if (n > 1) {
        size_t i;
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
    }
}


int compare_ints(const void * a, const void * b)
{
	const int *ia = (const int *) a;
	const int *ib = (const int *) b;
	
	return (*ia > *ib) - (*ia < *ib);
}


int compare_doubles (const void *a, const void *b){
	const double *da = (const double *) a;
	const double *db = (const double *) b;
 
	return (*da > *db) - (*da < *db);
}

/*=== Compare c[a] against c[b], where a and b are just indices. */
int compare_idxs_into_double_array (const void *a, const void *b){
	const int *ia = (const int *) a;
	const int *ib = (const int *) b;
	
	return (ptr_to_double_array_for_qsort[*ia] > ptr_to_double_array_for_qsort[*ib]) - (ptr_to_double_array_for_qsort[*ia] < ptr_to_double_array_for_qsort[*ib]);
}


void printMatrixDouble(double* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%lf ",X[i+nRows*j]);}
        printf(">\n");}
}

void printMatrixInt(int* X,int nRows, int nCols)
{
    int i,j;
    
    for(i = 0; i < nRows; i++) {
        printf("< ");
        for(j = 0; j < nCols; j++) {
            printf("%d ",X[i+nRows*j]);}
        printf(">\n");}
}



/* Kaplan-Meier estimator. If the highest value is censored K-L is
 * undefined for values above that.
 * For computing the mean we need an estimate, though. Fortunately, we
 * have an upper bound, and we assume that there are no deaths until
 * that upper bound.
 */
double kaplan_meier_mean(const int N, const int M, const double* y_in, const double* c_in, const double upper_bound){
    int c_idx, N_i, *rows, i, num_rows, d_i;
    double prod, km_mean, km_median, t_i, *y, *c; /* t_i instead of t[i] since we don't return the times*/
    bool median_set;
    if(N+M==0) {
        mexErrMsgTxt("Kaplan-Meier estimator undefined for empty population.");
    }

    /* Allocate Memory for Auxiliary Arrays. */
    rows = mxCalloc(N,sizeof(int));
    y = mxCalloc(N,sizeof(double));
    c = mxCalloc(M,sizeof(double));

    for(i=0; i<N; i++) y[i] = y_in[i];
    for(i=0; i<M; i++) c[i] = c_in[i];
    d_quick_sort(y, N);
    d_quick_sort(c, M);

    c_idx = 0;
    prod = 1.0;
    km_mean = 0;
/*    t=[];
    F=[];*/
    median_set = false;
    N_i = N+M;

    if (N>0){
        num_rows = 0;
        for(i=0; i<N-1; i++){
            if(y[i] + 1e-10 < y[i+1]){
                rows[num_rows++] = i;
            }
        }
        rows[num_rows++] = N-1;

        for(i=0; i<num_rows; i++){
            t_i = y[rows[i]];
            while (c_idx < M && c[c_idx] < t_i){
                N_i--;
                c_idx++;
            }
            if (i==0){
                d_i = rows[i]+1;
                km_mean += prod * y[rows[i]];
            } else {
                d_i = rows[i] - rows[i-1];
                km_mean += prod * (y[rows[i]] - y[rows[i-1]]);
            }
            prod *= (N_i-d_i+0.0)/(N_i+0.0);
/*            printf("N_i=%d, d_i=%d, after: y[rows[i]]=%lf, prod=%lf\n",N_i, d_i, y[rows[i]],prod); */
            N_i = N_i-d_i;
            /*F(i) = prod;*/
            if (!median_set && prod < 0.5 + 1e-10){
                if (prod < 0.5 - 1e-10){
                    km_median = t_i;
                } else {
                    if (i<num_rows-1){
                        km_median = 0.5*t_i + 0.5*y[rows[i+1]];
                    } else {
                        km_median = 0.5*t_i + 0.5*upper_bound;
                    }
                }
                median_set = true;
            }
        }
        /* Deal with remaining censored values (if the highest value is
         * uncensored, then prod=0, so nothing happens then)
         */
        /*km_mean += prod * (upper_bound-y[N-1]);*/
        km_mean += prod * (upper_bound-y[N-1])/2;
    } else {
        km_mean = upper_bound;
    }
    
    if (!median_set) km_median = upper_bound;
	/*=== Free Memory for Auxiliary Arrays. */
	mxFree(rows);
    mxFree(y);
    mxFree(c);

    return km_mean;
}




double logrank_statistic(const int num_periods, const int* N, const int* O, const int* N_1, const int* O_1, double* E, double* V){
    /* Compute logrank statistic for the specified "numbers at risk" (N),
     * observed events (O), "numbers at risk" in group 1 (N_1), and observed
     * events in group 1 (O_1)
     *
     * If N(i) <= 1, we don't include it in the sum (undefined variance, 0/0)
     */

    int p;
    double sum_V=0, denominator, numerator=0;

    for (p=0; p<num_periods; p++){
        E[p] = O[p] * ((N_1[p]+0.0)/N[p]);
    }

    for (p=0; p<num_periods; p++){
        if( N[p]<= 1 ){
            V[p] = 0;        
        } else {
            V[p] = E[p] * (1 - ((N_1[p]+0.0)/N[p])) * (N[p]-O[p]) / (N[p]-1.0);
            sum_V += V[p];
            numerator += (O_1[p]-E[p]);
        }
    }

    denominator = sqrt(sum_V);
    if (fabs(denominator) < 1e-6){
        if (fabs(numerator) < 1e-6){
/*            printf("result = 0\n");*/
            return 0;
        } else {
            
            printf("logrank num_periods=%d\nN:\n", num_periods);
            printMatrixInt(N, num_periods, 1);  

            printf("O:\n");
            printMatrixInt(O, num_periods, 1);  

            printf("N_1:\n");
            printMatrixInt(N_1, num_periods, 1);  

            printf("O_1:\n");
            printMatrixInt(O_1, num_periods, 1);  

            printf("E:\n");
            printMatrixDouble(E, num_periods, 1);  

            printf("V:\n");
            printMatrixDouble(V, num_periods, 1);  
            
            printf("numerator = %lf, denominator=%lf\n", numerator, denominator);
            mexErrMsgTxt("Division by zero in function logrank_statistic.");
        }
    }
    
/*    printf("result = %lf\n", numerator/denominator);*/
    return numerator/denominator;
}

                         
void Rcritval_cat_logrank(const double *x, const int nX, const double* y, const int* cens, const int *rows, const int nRows, const int *period, const int num_periods, const int *N_all, const int *O_all, const double kappa, double *critval_res, int* xleft, int* xright, int* numLeftPointer, int* numRightPointer){
	/* Declare Variables */
    int i, j, group_no, low, high, len, maxperiod, p, N_scalar, maxnumlocs, maxloc, s, numNonCens, numCens;
    int *N_1, *O_1, *N_add, *idx, *sorder, *maxlocs;
    double *logrank_stat, *km_mean, *allx, *E, *V, *y_in, *c_in;
    
	/*=== Allocate Memory for Auxiliary Arrays. */
    N_1 = mxCalloc(num_periods, sizeof(int));
    O_1 = mxCalloc(num_periods, sizeof(int));
    N_add = mxCalloc(num_periods, sizeof(int));
    idx = mxCalloc(nX, sizeof(int));
    sorder = mxCalloc(nRows+1, sizeof(int));
    maxlocs = mxCalloc(nRows,sizeof(int));
    
    logrank_stat = mxCalloc(nRows, sizeof(double));
    km_mean = mxCalloc(nRows+1,sizeof(double));
    allx = mxCalloc(nRows+1,sizeof(double));
    E = mxCalloc(num_periods, sizeof(double));
    V = mxCalloc(num_periods, sizeof(double));
    y_in = mxCalloc(nX, sizeof(double));
    c_in = mxCalloc(nX, sizeof(double));
/*    
    printf("nRows=%d, rows:\n", nRows);
	printMatrixInt(rows, 1, nRows);  
*/
    
    /* Get Kaplan-Meier estimates of the mean & median of each group. */
    for(i=0; i<=nRows; i++){
        if(i==0){
            low = 0;
        } else {
            low = rows[i-1];
        }
        
        if(i==nRows){
            high = nX;
        } else {
            high = rows[i];
        }

        numNonCens = 0;
        numCens = 0;
        for(j=low; j<high; j++){
            if(cens[j]==0){
                y_in[numNonCens++] = y[j];
            } else {
                c_in[numCens++] = y[j];
            }
        }
        km_mean[i] = kaplan_meier_mean(numNonCens, numCens, y_in, c_in, kappa);
/*
        if (i==1){
            printf("Computing km_mean[%d], numNonCens=%d, numCens=%d, kappa=%lf, y_in:\n", i, numNonCens, numCens, kappa);
            printMatrixDouble(y_in, numNonCens, 1);  
            printf("c_in:\n");
            printMatrixDouble(c_in, numCens, 1);  
            printf("result = %lf:\n", km_mean[i]);
        }
*/
    }

	/*=== Sort categories by km_mean. */
    for(i=0; i<nRows+1; i++){
        sorder[i] = i;
    }
	ptr_to_double_array_for_qsort = km_mean;
    dp_quick_sort(sorder, nRows+1);
    
    for (i=0; i<num_periods; i++){
        N_1[i] = 0;
        O_1[i] = 0;
    }
/*
    printf("km_mean:\n");
   	printMatrixDouble(km_mean, nRows+1, 1);  
    
    printf("sorder:\n");
   	printMatrixInt(sorder, nRows+1, 1);  
*/
    
/*
    printf("nRows=%d, rows:\n", nRows);
    printMatrixInt(rows, 1, nRows);  
*/    
    /*=== For adding one group at a time, compute resulting logrank statistic. */
    for (s=0; s<nRows; s++){ /* for debug: group_no<nRows+1 */
        group_no = sorder[s];
    
        /* Get indices of group_no'th group. */
        if (group_no == 0){
            low = 0;
        } else {
            low = rows[group_no-1];
        }
        if (group_no == nRows){
            high = nX;
        } else {
            high = rows[group_no];;
        }
        len = high-low;
        
        for(i=0; i<len; i++){
            idx[i] = low+i;
        }
        
        maxperiod = -1;
        for (i=0; i<len; i++){
            maxperiod = MAX(maxperiod, period[idx[i]]);
        }
        
        for (p=0; p<num_periods; p++){
            N_add[p] = 0;
        }
        
        /* Collect observations in O_1, and helper N_add to collect numbers at risk, N_1 */
        for (i=0; i<len; i++){
            p = period[idx[i]]-1;
			if( p>=0 ){
				if(cens[idx[i]]==0){
					O_1[p]++;
				}
				N_add[p]++;
			}
        }
/*
        printf("N_add:\n");
        printMatrixInt(N_add, 1, num_periods);  
*/        
        
        N_scalar=0;
        for(p=num_periods-1; p>=0; p--){
            N_scalar += N_add[p];
            N_1[p] += N_scalar;
            if(N_1[p] < O_1[p]){
                printf("x\n");
                printMatrixDouble(x, 1, nX);  
                
                printf("rows\n");
                printMatrixInt(rows, 1, nRows);  
                
                printf("y\n");
                printMatrixDouble(y, 1, nX);  
                
                printf("period:\n");
                printMatrixInt(period, 1, nX);  

                printf("cens:\n");
                printMatrixInt(cens, 1, nX);  
                
                printf("N_all:\n");
                printMatrixInt(N_all, 1, num_periods);  

                printf("O_all:\n");
                printMatrixInt(O_all, 1, num_periods);  

                printf("km_mean:\n");
                printMatrixDouble(km_mean, 1, nRows+1);  
                
                printf("s=%d, group_no=%d, sorder:\n", s, group_no);
                printMatrixInt(sorder, 1, nRows+1);  
                
                printf("idx, len=%d:\n", len);
                printMatrixInt(idx, 1, len);  
                
                printf("x[idx]:");
                for(i=0; i<len; i++){
                    printf(" %lf ", x[idx[i]]);
                }
                printf("\n");

                printf("cens[idx]:");
                for(i=0; i<len; i++){
                    printf(" %d ", cens[idx[i]]);
                }
                printf("\n");

                printf("maxperiod = %d, N_add:\n", maxperiod);
                printMatrixInt(N_add, 1, num_periods);  
                
                printf("N_1:\n");
                printMatrixInt(N_1, 1, num_periods);  

                printf("O_1:\n");
                printMatrixInt(O_1, 1, num_periods);  
                
                mexErrMsgTxt("O_1[p] > N_1[p] -- should be impossible");
            }
        }

/*        printf("maxperiod = %d, N_scalar=%d\n", maxperiod, N_scalar);
        if( group_no==500 ){
            printf("logrank num_periods=%d\nN:\n", num_periods);
            printMatrixInt(N_all, 1, num_periods);  

            printf("O:\n");
            printMatrixInt(O_all, 1, num_periods);  

            printf("N_1:\n");
            printMatrixInt(N_1, 1, num_periods);  

            printf("O_1:\n");
            printMatrixInt(O_1, 1, num_periods);  
            
            return;
        }
 */       
        logrank_stat[s] = logrank_statistic(num_periods, N_all, O_all, N_1, O_1, E, V);
    }
    
    critval_res[0] = -1;
    for(i=0; i<nRows; i++){
        if(fabs(logrank_stat[i]) > critval_res[0]-1e-6){
            if(fabs(logrank_stat[i]) > critval_res[0]+1e-6){
                critval_res[0] = fabs(logrank_stat[i]);
                maxnumlocs = 0;
            }
            maxlocs[maxnumlocs] = i;
            maxnumlocs = maxnumlocs + 1;
        }
    }
/*
	printf("\nlogrank_stat:\n");
	printMatrixDouble(logrank_stat, nRows, 1);  
	printf("\n\n");

	printf("\nmaxlocs:\n");
	printMatrixInt(maxlocs, maxnumlocs, 1);  
	printf("\n\n");
*/
	maxloc = maxlocs[rand()%maxnumlocs];
/*	printf("%maxnumlocs = %d\n", maxnumlocs);
//	printf("maxloc %d\n", maxloc); */

/*	maxloc = maxlocs[0]; /* TODO: remove this after debugging */
	
    /* Take one x value from each unique set */
    for(i=0; i<nRows; i++){
        allx[i] = x[rows[i]-1];
    }
    allx[nRows] = x[nX-1];
    
    numLeftPointer[0] = maxloc+1;
	numRightPointer[0] = nRows+1-maxloc-1;

/* 
    printf("%allx:\n");
   	printMatrixDouble(allx, nRows+1, 1);  
*/
    
	/* Now we can fill the result arrays xleft and xright as usual. */
	for(i=0; i<maxloc+1; i++){
	    xleft[i] = (int) allx[sorder[i]];
	}
	for(i=maxloc+1; i<nRows+1; i++){
	    xright[i-maxloc-1] = (int) allx[sorder[i]];
	}

/*
	printf("\nxleft  inside:\n");
	printMatrixInt(xleft, maxloc+1, 1);  
	printf("\n\n");

	printf("\nxright inside:\n");
	printMatrixInt(xright, nRows+1-maxloc-1, 1);  
	printf("\n\n");
*/

	/*=== Sort outputs. */
	quick_sort(xleft, maxloc+1);
	quick_sort(xright, nRows+1-maxloc-1);

    
   	/*=== Free Memory for Auxiliary Arrays. */
    mxFree(N_1);
    mxFree(O_1);
    mxFree(N_add);
    mxFree(idx);
    mxFree(maxlocs);
    mxFree(logrank_stat);
    mxFree(km_mean);
    mxFree(allx);
    mxFree(E);
    mxFree(V);
}

void Rcritval_cont_logrank(const double *x, const int nX, const int* cens, const int *rows, const int nRows, const int *period, const int num_periods, const int *N_all, const int *O_all, double *critval_res, double *cutval_res){
	/* Declare Variables */
    int i, group_no, low, high, len, maxperiod, p, N_scalar, maxnumlocs, maxloc, cutloc;
    int *N_1, *O_1, *N_add, *idx, *maxlocs;
    double u,*logrank_stat, *E, *V;
    
	/*=== Allocate Memory for Auxiliary Arrays. */
    N_1 = mxCalloc(num_periods, sizeof(int));
    O_1 = mxCalloc(num_periods, sizeof(int));
    N_add = mxCalloc(num_periods, sizeof(int));
    idx = mxCalloc(nX, sizeof(int));
    maxlocs = mxCalloc(nRows,sizeof(int));
    
    logrank_stat = mxCalloc(nRows, sizeof(double));
    E = mxCalloc(num_periods, sizeof(double));
    V = mxCalloc(num_periods, sizeof(double));
/*    
    printf("nRows=%d, rows:\n", nRows);
	printMatrixInt(rows, 1, nRows);  
*/
    for (i=0; i<num_periods; i++){
        N_1[i] = 0;
        O_1[i] = 0;
    }
/*
    printf("nRows=%d, rows:\n", nRows);
    printMatrixInt(rows, 1, nRows);  
*/    
    /*=== For adding one group at a time, compute resulting logrank statistic. */
    for (group_no=0; group_no<nRows; group_no++){ /* for debug: group_no<nRows+1 */
        /* Get indices of group_no'th group. */
        if (group_no == 0){
            low = 0;
        } else {
            low = rows[group_no-1];
        }
        if (group_no == nRows){
            high = nX; /* only during debugging, when group_no gets that high */
        } else {
            high = rows[group_no];;
        }
        len = high-low;
        
        for(i=0; i<len; i++){
            idx[i] = low+i;
        }        

        maxperiod = -1;
        for (i=0; i<len; i++){
            maxperiod = MAX(maxperiod, period[idx[i]]);
        }

        for (p=0; p<maxperiod; p++){
            N_add[p] = 0;
        }
        
        /* Collect observations in O_1, and helper N_add to collect numbers at risk, N_1 */
        for (i=0; i<len; i++){
            p = period[idx[i]]-1;
			if( p>=0 ){
				if(cens[idx[i]]==0){
					O_1[p]++;
				}
				N_add[p]++;
			}
        }
/*
        printf("N_add:\n");
        printMatrixInt(N_add, 1, num_periods);  
*/        
        
        N_scalar=0;
        for(p=maxperiod-1; p>=0; p--){
            N_scalar += N_add[p];
            N_1[p] += N_scalar;
        }

/*        printf("maxperiod = %d, N_scalar=%d\n", maxperiod, N_scalar);
        if( group_no==500 ){
            printf("logrank num_periods=%d\nN:\n", num_periods);
            printMatrixInt(N_all, 1, num_periods);  

            printf("O:\n");
            printMatrixInt(O_all, 1, num_periods);  

            printf("N_1:\n");
            printMatrixInt(N_1, 1, num_periods);  

            printf("O_1:\n");
            printMatrixInt(O_1, 1, num_periods);  
            
            return;
        }
 */       
        logrank_stat[group_no] = logrank_statistic(num_periods, N_all, O_all, N_1, O_1, E, V);
    }
    
    critval_res[0] = -1;
    for(i=0; i<nRows; i++){
        if(fabs(logrank_stat[i]) > critval_res[0]-1e-10){
            if(fabs(logrank_stat[i]) > critval_res[0]+1e-10){
                critval_res[0] = fabs(logrank_stat[i]);
                maxnumlocs = 0;
            }
            maxlocs[maxnumlocs] = i;
            maxnumlocs = maxnumlocs + 1;
        }
    }
/*
	printf("\nlogrank_stat:\n");
	printMatrixDouble(logrank_stat, nRows, 1);  
	printf("\n\n");

	printf("\nmaxlocs:\n");
	printMatrixInt(maxlocs, maxnumlocs, 1);  
	printf("\n\n");
*/
	maxloc = maxlocs[rand()%maxnumlocs];
/*	printf("%maxnumlocs = %d\n", maxnumlocs);
//	printf("maxloc %d\n", maxloc); */

/*	maxloc = maxlocs[0]; // TODO: remove this after debugging */
	
	/*=== Get cutval. */
	cutloc = rows[maxloc]-1;
    
    u = rand()/(RAND_MAX+0.0);
/*	printf("u = %lf\n", u);*/
/*	cutval_res[0] = (x[cutloc] + x[cutloc+1])/2;*/
	cutval_res[0] = ((1-u)*(x[cutloc]+1e-4) + u*(x[cutloc+1]-1e-4));

   	/*=== Free Memory for Auxiliary Arrays. */
    mxFree(N_1);
    mxFree(O_1);
    mxFree(N_add);
    mxFree(idx);
    mxFree(maxlocs);
    mxFree(logrank_stat);
    mxFree(E);
    mxFree(V);
}


void  buildTheTree(const double* X, const double* y, const int* cens, const int SplitMinUncens, const int numFeaturesType1, const double regtree_p, const double percentageFeatures, const int* iscat, const int N, const int nvars, const int log10transformed, const double kappa, const double cutoff_penalty_factor, int* nodenumber, int* parent, mxArray* ysub, mxArray* censsub, int* cutvar, double* cutpoint, int* leftchildren, int* rightchildren, int* nodesize, mxArray* catsplit, int* numNodesPointer, int* numNcatsplitPointer){
	int i, j, k, offset, nRest, ncatsplit, tnode, nextunusednode, Nnode, bestvar, nRandom, numVarsToConsider, jvar, numrows, xcat, numBestLeft, numBestRight, nleft, nright, currnode, parent_node, catsplit_index, num_compatible, num_missing_to_left, num_missing_to_right, numUncensNode;
	int *censnode, *censnode_idx, *noderows, *leftside, *rightside, *assignednode, *randomPermutation, *idx, *rows;
    int *period, *period_idx, *N_all, *O_all, *C_all;
	int *xleft, *xright, *numLeftPointer, *numRightPointer, *bestleft, *bestright, *xleftForResult, *xrightForResult;
    int num_periods, curr_period, p;
	double ybar, sst, mincost, bestcrit, bestcut, probForceSplitOnFeatureType1;
	double *xnoderow, *ynode, *ynode_idx, *x, *critvalPointer, *cutvalPointer, *ysub_for_result;
    int *censsub_for_result;
    double *period_time;
    double last_time, this_time;
	bool impure, ismember;
	int dims_left[2], dims_right[2], dims[2];
	mxArray *mx_xleft, *mx_xright, *mx_ysub, *mx_censsub;
    
	/*=== Allocate Memory for Auxiliary Arrays.*/
    noderows = mxCalloc(N,sizeof(int));
    leftside = mxCalloc(N,sizeof(int));
    rightside = mxCalloc(N,sizeof(int));
    assignednode = mxCalloc(N,sizeof(int));
    idx = mxCalloc(N,sizeof(int));
    rows = mxCalloc(N,sizeof(int));
    bestleft = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
    bestright = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
    censnode = mxCalloc(N,sizeof(int));
    censnode_idx = mxCalloc(N,sizeof(int));
    period = mxCalloc(N,sizeof(int));
    period_idx = mxCalloc(N,sizeof(int));
    N_all = mxCalloc(N,sizeof(int));
    C_all = mxCalloc(N,sizeof(int));
    O_all = mxCalloc(N,sizeof(int));

    randomPermutation = mxCalloc(nvars,sizeof(int));
    
	xnoderow = mxCalloc(N,sizeof(double));
	ynode = mxCalloc(N,sizeof(double));
    ynode_idx = mxCalloc(N,sizeof(double));
    x = mxCalloc(N,sizeof(double));
    period_time = mxCalloc(N+1,sizeof(double));
    
	/*=== Allocate Memory for Auxiliary Arrays for calling Rcritval_cat & Rcritval_cont.*/
	xleft = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
	xright = mxCalloc(N,sizeof(int)); /* possible values for a var limited by the total # training data */
	numLeftPointer = mxCalloc(1,sizeof(int));
	numRightPointer = mxCalloc(1,sizeof(int));
	critvalPointer = mxCalloc(1,sizeof(double));
	cutvalPointer = mxCalloc(1,sizeof(double));

/*
	printf("X:\n");
	printMatrixDouble(X, N, nvars);
	printf("\n\n");
*/	
	/*=== Initialize variables. */
	ncatsplit = 0;
	nextunusednode = 2;
	nodenumber[0] = 1;
	for(i=0; i<N; i++){
		assignednode[i] = 1;
	}
    
	/*=== Keep processing nodes until done. */
	tnode = 1;
/* 	while(tnode < 2){ for debug */
	while(tnode < nextunusednode){
        /*=== Record information about this node.*/
		/* Matlab: noderows = find(assignednode==tnode); */
		Nnode = 0;
		for(i=0; i<N; i++){
	       if( assignednode[i] == tnode ){
	           noderows[Nnode] = i;
	           Nnode = Nnode + 1;
	       }
		}
        /*printf("tnode = %d\n", tnode);*/

		/* Matlab: ynode = y(noderows); */
        numUncensNode = 0;
		for(i=0; i<Nnode; i++){
			ynode[i] = y[noderows[i]];
            censnode[i] = cens[noderows[i]];
            numUncensNode += (1-censnode[i]);
/*            printf("ynode[%d]=%lf\n",i,ynode[i]);*/
		}
/*		
		printf("censnode:\n");
		printMatrixInt(censnode, 1, Nnode);
		printf("\n\n");
        printf("tnode = %d, numUncensNode=%d\n", tnode, numUncensNode);
  
		printf("ynode:\n");
		printMatrixDouble(ynode, 1, Nnode);
		printf("\n\n");
*/		
		/*=== Compute mean, variance and related statistics for this node. */
		/* Matlab: ybar = mean(ynode); */
		ybar = 0;
		for(i=0; i<Nnode; i++){
       		ybar += ynode[i]; 
		}
		ybar = ybar/Nnode;

		/* Matlab: sst = norm(ynode-ybar)^2;   % total sum of squares at this node */
		sst = 0;
		for(i=0; i<Nnode; i++){
			sst = sst + (ynode[i]-ybar)*(ynode[i]-ybar);
		}

        if (Nnode > 1) {
            mincost = sqrt(sst / (Nnode-1)); /* stddev of ynode */
        } else {
            mincost = 0;
        }
		impure = (mincost > 0);

/*		printf("sst=%lf, mincost = %lf, impure=%d \n\n", sst, mincost, impure?1:0);*/
		
		/*=== Initialize variables before looping over possible split vars. */
		bestcrit          = -1e12;
		nodesize[tnode-1]   = Nnode;
		cutvar[tnode-1]     = 0;
		cutpoint[tnode-1]   = 0;
		leftchildren[tnode-1] = 0;
		rightchildren[tnode-1] = 0;

        /*=== Consider splitting this node. */
		if ( (numUncensNode >= SplitMinUncens) && impure ){     /* split only impure nodes with more than a threshold of uncensored values*/
			/* Matlab: Xnode = X(noderows,:);  I don't want to deal with temporary matrices, so I'll work around by indexing. */
			
			bestvar = -1;
			bestcut = 0;

			/*=== First decision: force split on algorithm parameter? */
            probForceSplitOnFeatureType1 = pow( (log2(N)-log2(Nnode)) / (log2(N)-1), regtree_p);
/*            probForceSplitOnFeatureType1 = 0.5; /*MAX(probForceSplitOnFeatureType1, 0.5);*/
            /*probForceSplitOnFeatureType1 = 0.8;*/
			/*      //probForceSplitOnFeatureType1= (p ^(1/(N-2.0)))^(Nnode-2); */

/*            
            printf("Nnode=%d, N=%d, probForceSplitOnFeatureType1=%lf\n",Nnode, N, probForceSplitOnFeatureType1);
*/            
			if (probForceSplitOnFeatureType1*RAND_MAX > rand()){
				/*=== Force split on algorithm parameter. */
				nRandom = numFeaturesType1;
                nRest = nvars-numFeaturesType1;
 /*               printf("forceSplitOnFeatureType1 on\n");*/
			} else {
				/*=== Don't force any type of split variable. */
				nRandom = nvars;
                nRest=0;
/*                printf("forceSplitOnFeatureType1 off\n");*/
			}
			
			/* Matlab: randomPermutation = randperm(nRandom); */
			for(i=0; i<nvars; i++){
				randomPermutation[i] = i;
			}
			shuffle(randomPermutation, nRandom);    /* take out for debugging */ 

            if (nRest>0){
                shuffle(randomPermutation+nRandom, nRest); 
/*                printf("randperm forced:\n");*/
            } else {
/*                printf("randperm not forced:\n");*/
            }
/*            printMatrixInt(randomPermutation, nvars, 1); */


/*			printf("percentageFeatures=%lf\n",percentageFeatures);*/
			numVarsToConsider = MAX(1, (int) floor(percentageFeatures*nRandom));
/*            printf("nvars=%d, numVarsToConsider=%d\n",nvars,numVarsToConsider);*/


            ptr_to_double_array_for_qsort = ynode;
            for(j=0; j<Nnode; j++){
                idx[j] = j;
            }
            dp_quick_sort(idx, Nnode);
                        
            /* First, compute period for uncensored data. (Censored data with same y value could be before or after the first uncensored one in the ordering, thus we need 2 passes.) */
            num_periods=0;
            last_time = -1e10;
            for(j=0; j<Nnode; j++){
                i = idx[j];
/*                printf("j=%d, i=%d\n", j, i);*/
                if( censnode[i] == 0 ){
                    if( ynode[i] > last_time + 1e-6 ){
                        num_periods++;
                        period_time[num_periods] = ynode[i];
                        last_time = ynode[i];
                    }
                    period[i] = num_periods; /* Matlab indexing, period starts at 1 */
                }
/*                printf("period[i]=%d\n", period[i]);*/
            }
/*
            printf("num_periods=%d, period:\n", num_periods);
            printMatrixInt(period, Nnode, 1);
            printf("\n\n"); 
*/            
            /* Then, fill in period for censored data. */
            curr_period = 0;
            for(i=0; i<Nnode; i++){
                if( censnode[i]==1 ){
                    if( curr_period == num_periods ){
                        this_time = 1e10;
                    } else {
                        this_time = period_time[curr_period+1];
                    }
                    while ( ynode[i] > this_time-1e-6 ){ /* i.e. >= */
                        curr_period++;
                        if( curr_period == num_periods ){
                            this_time = 1e10;
                        } else {
                            this_time = period_time[curr_period+1];
                        }
                    }
                    period[i] = curr_period;
                }
            }
/*            
            printf("num_periods=%d, period:\n", num_periods);
            printMatrixInt(period, Nnode, 1);
            printf("\n\n"); 
*/
            /* Compute initial logrank counters when all data in 2nd group. */
            for(p=0; p<num_periods; p++){
                O_all[p] = 0;
                C_all[p] = 0;
            }

            /* Go through data points, and collect the events. */
            for(i=0; i<Nnode; i++){
                p = period[i];
                if(censnode[i] == 0){
                    O_all[p-1]++;
                } else {
                    if (p>0){
                        C_all[p-1]++;
                    } /* else the point counts as if it hadn't even happened */
                }
            }

            /* Go through periods, and collect N. */
            N_all[0] = Nnode;
            for(p=1; p<num_periods; p++){
                N_all[p] = N_all[p-1] - O_all[p-1] - C_all[p-1];
            }
            /* End of this part of special logrank code. */
/*
            printf("period, Nnode=%d:\n", Nnode);
            printMatrixInt(period, 1, Nnode);
            printf("\n\n"); 
            
            printf("N_all, num_periods=%d:\n", num_periods);
            printMatrixInt(N_all, num_periods, 1);
            printf("\n\n"); 
            
            printf("O_all, num_periods=%d:\n", num_periods);
            printMatrixInt(O_all, num_periods, 1);
            printf("\n\n"); 
*/
			/*=== Second decision: which one? */
			for(i=0; i<nvars; i++){ /* we allow numVarsToConsider options, but must have at least one real option */
				jvar=randomPermutation[i]; 
/*                jvar = i; /* TODO: remove this after debugging. */
				xcat = iscat[jvar];
			
				/* Matlab: Xnoderow = Xnode(:,jvar) */
                offset = jvar*N; /* index into matrix: row + column*numRows */
				for(j=0; j<Nnode; j++){
					xnoderow[j] = X[noderows[j] + offset];  
				}
/*
                printf("jvar=%d, xnoderow:\n", jvar);
				printMatrixDouble(xnoderow, Nnode, 1);
				printf("\n\n"); 
*/
                /* Matlab: xnoderow=xnoderow(1:Nnode); %only sort first Nnode elements.
				// Matlab: [x,idx] = sort(xnoderow);          % get sorted jth x variable */
				ptr_to_double_array_for_qsort = xnoderow;
			    for(j=0; j<Nnode; j++){
			        idx[j] = j;
			    }
                dp_quick_sort(idx, Nnode);
                for(j=0; j<Nnode; j++){
			        x[j] = xnoderow[idx[j]];
			    }
/*
                printf("jvar=%d, x:\n", jvar);
				printMatrixDouble(x, Nnode, 1);
				printf("\n\n");
*/
                /*=== Determine if there's anything to split along this variable. */
				if (x[Nnode-1]-x[0] < 1e-10){
					continue;
				}
/*
            printf("jvar=%d, x\n", jvar);
            printMatrixDouble(x, 1, Nnode);
            printf("\n\n"); 
*/				
				/* Matlab: rows = find(x(1:end-1)+maxeps < x(2:end));
				// WATCH OUT: rows holds the indices as original in Matlab, not C style (but it itself is referenced standard C style starting at 0) */
				numrows = 0;
				for (j=0; j<Nnode-1; j++){
					if (x[j]+1e-10 < x[j+1]){ 
						rows[numrows++] = j+1; /* the +1 here is to make this compatible with calling Rcritval_cat & cont from Matlab directly. */
					}
				}
				if (numrows==0){
					continue;
				}
/*
            printf("jvar=%d, cat=%d, rows, numrows=%d:\n", jvar, xcat, numrows);
            printMatrixInt(rows, 1, numrows);
*/                
                for(j=0; j<Nnode; j++){
                    ynode_idx[j] = ynode[idx[j]];
                    censnode_idx[j] = censnode[idx[j]];
                    period_idx[j] = period[idx[j]];
                }
            
				/*=== Do the core work: get the best split of the variable and its quality. */
				if (xcat>0){
					/* Rcritval_cat(x, ycum, rows, Nnode, numrows, critvalPointer, xleft, xright, numLeftPointer, numRightPointer); */
                    Rcritval_cat_logrank(x, Nnode, ynode_idx, censnode_idx, rows, numrows, period_idx, num_periods, N_all, O_all, kappa, critvalPointer, xleft, xright, numLeftPointer, numRightPointer);
/*                    printf("cat var %d, critval=%lf\n",jvar,critvalPointer[0]);*/

/*
					printf("numleft=%d, xleft:\n", numLeftPointer[0]);
					printMatrixInt(xleft, numLeftPointer[0], 1);
					printf("\n\n");

					printf("numright=%d, xright:\n", numRightPointer[0]);
					printMatrixInt(xright, numRightPointer[0], 1);
					printf("\n\n");
*/
				} else {
					/* Rcritval_cont(x, ycum, rows, Nnode, numrows, critvalPointer, cutvalPointer); */
                    Rcritval_cont_logrank(x, Nnode, censnode_idx, rows, numrows, period_idx, num_periods, N_all, O_all, critvalPointer, cutvalPointer);                    
/*                    printf("cont var %d, critval=%lf\n",jvar,critvalPointer[0]);*/
				}
/*            printf("\n\n"); */

				/*=== Change best split if this one is best so far. */
				if (critvalPointer[0] > bestcrit + 1e-10){
					bestcrit = critvalPointer[0];
					bestvar = jvar;
					if (xcat>0){
						numBestLeft = numLeftPointer[0];
						numBestRight = numRightPointer[0];
						for(j=0; j<numBestLeft; j++){
							bestleft[j] = xleft[j];
						}
						for(j=0; j<numBestRight; j++){
							bestright[j] = xright[j];
						}
					} else {
						bestcut = cutvalPointer[0];
					}
				}
/*
                printf("i=%d, jvar = %d, xcat = %d, bestcrit = %lf, bestvar = %d, critval = %lf\n", i, jvar, xcat, bestcrit, bestvar, critvalPointer[0]);
*/
/*
				printf("jvar=%d, bestvar=%d, bestcrit=%lf:\n", jvar, bestvar, bestcrit);
				if (iscat[bestvar]){
					printf("numBestLeft=%d, bestleft:\n", numBestLeft);
					printMatrixInt(bestleft, numBestLeft, 1);
					printf("\n\n");

					printf("numBestRight=%d, bestright:\n", numBestRight);
					printMatrixInt(bestright, numBestRight, 1);
					printf("\n\n");
				} else {
					printf("bestcut=%lf\n", bestcut);
				}
*/				
				if (i >= numVarsToConsider-1 && bestcrit > -1e11){
/*                    printf("Have looked at %d variables, bestcrit=%lf, breaking\n",i+1,bestcrit);*/
					break;
				}
/*                printf("i=%d, bestcrit=%lf\n",i,bestcrit);*/
			}
/*
            printf("\n");
*/			
			/*=== Split this node using the best rule found. */
			if (bestvar == -1){
				/* Terminal node
		        // printf("Terminal node with %d data points and impure=%d\n", Nnode, impure?1:0); */
        	} else {
				for (j=0; j<Nnode; j++){
					x[j] = X[noderows[j] + bestvar*N];
				}
/*				
				printf("\nnoderows:\n");
				printMatrixInt(noderows, Nnode, 1);  
				printf("\n\n");
				
				printf("\nx to be split:\n");
				printMatrixDouble(x, Nnode, 1);  
				printf("\n\n");
*/				
    
				if (iscat[bestvar]){
/*				printf("splitting on cat %d\n", bestvar); */
					cutvar[tnode-1] = -(bestvar+1);          /* negative indicates cat. var. split */
					ncatsplit++;  	   /* index into catsplit cell array */
					cutpoint[tnode-1] = ncatsplit;

                    /* 1: To get all compatible values, walk up the tree, looking
                     * for a split on the same parameter. If none is found
                     * take the initial domain of that parameter. */
                    currnode = tnode;
                    while (currnode > 1){
/*                        printf("currnode = %d\n", currnode);*/
                        parent_node = parent[currnode-1];
/*                        printf("parent_node = %d, cutvar[parent_node-1]=%d\n", parent_node, cutvar[parent_node-1]);*/
                        if (cutvar[parent_node-1] == -(bestvar+1)){
                            /* Take values from there, depending on whether which child we are */
                            catsplit_index = cutpoint[parent_node-1];
/*                            printf("catsplit_index = %d\n", catsplit_index);*/

                            if (leftchildren[parent_node-1] == currnode){
                                mx_to_get_compatible_values = mxGetCell(catsplit, catsplit_index-1);
                            } else {
                                if (! (rightchildren[parent_node-1] == currnode)){
                                    mexErrMsgTxt("currnode must either be left or right child of its parent.");
                                }
                                mx_to_get_compatible_values = mxGetCell(catsplit, catsplit_index-1+N);
                            }
                            break;
                        }
                        currnode = parent_node;
                    }
/*                    printf("done loop\n");*/
                    if (currnode == 1){
                        /* Get compatible values from initial domain. */
/*                        printf("bestvar=%d\n",bestvar);*/
                        mx_to_get_compatible_values = mxGetCell(domains_cat, bestvar);
                    }
                    /* Get compatible values from mx_to_get_compatible_values. */
                    num_compatible = mxGetNumberOfElements(mx_to_get_compatible_values);
                    compatible_values = (int*) mxGetData(mx_to_get_compatible_values);
                    
/*                    printf("num_compatible=%d\n",num_compatible);*/
                    /* 2: For each compatible but missing value choose a side u.a.r. */
                    missing_values_for_left = mxCalloc(num_compatible,sizeof(int));
                    missing_values_for_right = mxCalloc(num_compatible,sizeof(int));
                    num_missing_to_left = 0;
                    num_missing_to_right = 0;
                    for (i=0; i<num_compatible; i++){
                        for (j=0; j<numBestLeft; j++){
                            if (compatible_values[i] == bestleft[j]) break;
                        }
                        if (j == numBestLeft){
                            for (j=0; j<numBestRight; j++){
                                if (compatible_values[i] == bestright[j]) break;
                            }
                            if (j == numBestRight){
                                /* Missing but compatible value: choose side u.a.r. */
                                if (rand()%2 == 0){
                                    missing_values_for_left[num_missing_to_left++] = compatible_values[i];
                                } else {
                                    missing_values_for_right[num_missing_to_right++] = compatible_values[i];
                                }
                            }
                        }
                    }
/*                    printf("num_missing_to_left=%d\n",num_missing_to_left);
                    printf("num_missing_to_right=%d\n",num_missing_to_right);*/
                    
                    /* 3: Merge the determined and the randomly assigned missing values */
                    for (i=num_missing_to_left; i<num_missing_to_left+numBestLeft; i++){
                        missing_values_for_left[i] = bestleft[i-num_missing_to_left];
                    }
                    quick_sort(missing_values_for_left, num_missing_to_left+numBestLeft);

                    for (i=num_missing_to_right; i<num_missing_to_right+numBestRight; i++){
                        missing_values_for_right[i] = bestright[i-num_missing_to_right];
                    }
                    quick_sort(missing_values_for_right, num_missing_to_right+numBestRight);

                    /* 4: Put that information into the cell array. */

		/*=== Set up the structures to fill the cell array output. */
					dims_left[0] = 1;
					dims_left[1] = num_missing_to_left+numBestLeft; /* numBestLeft; */
					
					dims_right[0] = 1;
					dims_right[1] = num_missing_to_right+numBestRight; /* numBestRight; */
					
					mx_xleft = mxCreateNumericArray(2, dims_left, mxINT32_CLASS, mxREAL);
					mx_xright = mxCreateNumericArray(2, dims_right, mxINT32_CLASS, mxREAL);
					
					mxSetCell(catsplit, ncatsplit-1, mx_xleft);
					mxSetCell(catsplit, ncatsplit-1+N, mx_xright);
			
					xleftForResult = (int*) mxGetData(mx_xleft);
					xrightForResult = (int*) mxGetData(mx_xright);

					/*=== Copy result from our temporary arrays to the ones associated with the output. */
		    /*
					for(i=0; i<numBestLeft; i++){
						xleftForResult[i] = bestleft[i];
					}
					for(i=0; i<numBestRight; i++){
						xrightForResult[i] = bestright[i];
					}
                    */
                    for(i=0; i<num_missing_to_left+numBestLeft; i++){
						xleftForResult[i] = missing_values_for_left[i];
					}
					for(i=0; i<num_missing_to_right+numBestRight; i++){
						xrightForResult[i] = missing_values_for_right[i];
					}
                    mxFree(missing_values_for_left);
                    mxFree(missing_values_for_right);
	
					/* Matlab: leftside = ismember(x,bestleftrightcell{1});
					// Matlab: rightside = ismember(x,bestleftrightcell{2}); */
					nleft = 0;
					nright = 0;
/*
					printf("numBestLeft=%d\n", numBestLeft);
					printf("\nbestleft:\n");
					printMatrixInt(bestleft, numBestLeft, 1);  
					printf("\n\n");

					printf("numBestRight=%d\n", numBestRight);
					printf("\nbestright:\n");
					printMatrixInt(bestright, numBestRight, 1);  
					printf("\n\n");

					printf("\nx:\n");
					printMatrixDouble(x, Nnode, 1);  
					printf("\n\n");
*/
					for (j=0; j<Nnode; j++){
                        
						ismember = false;
						for (k=0; k<numBestLeft; k++){
							if ( ((int) floor(x[j]+0.5)) == bestleft[k]) ismember = true;
						}
						
						if (ismember){
							leftside[nleft] = j;
							nleft = nleft+1;
						} else {
							rightside[nright] = j;
							nright = nright+1;
						}
					}
				} else {

/*					printf("splitting on cont %d\n", bestvar); */
					cutvar[tnode-1] = bestvar+1;
					cutpoint[tnode-1] = bestcut;
					
					/* Matlab: leftside = x<=bestcut; %logical   
					// Matlab: rightside = ~leftside; */
					nleft = 0;
					nright = 0;
					for (j=0; j<Nnode; j++){
						if (x[j] <= bestcut){
							leftside[nleft] = j;
							nleft = nleft+1;
						} else {
							rightside[nright] = j;
							nright = nright+1;
						}
					}
				}
/*				
				printf("\nleftside:\n");
				printMatrixInt(leftside, nleft, 1);  
				printf("\n\n");

				printf("\nrightside:\n");
				printMatrixInt(rightside, nright, 1);  
				printf("\n\n");
*/			

				leftchildren[tnode-1] = nextunusednode;
				rightchildren[tnode-1] = nextunusednode+1;
				for (j=0; j<nleft; j++){
					assignednode[noderows[leftside[j]]] = nextunusednode;
				}
				for (j=0; j<nright; j++){
					assignednode[noderows[rightside[j]]] = nextunusednode+1;
				}
	
/*				
				printf("\nassignednode:\n");
				printMatrixInt(assignednode, N, 1);  
				printf("\n\n");
*/				

				nodenumber[nextunusednode-1] = nextunusednode;
				nodenumber[nextunusednode-1+1] = nextunusednode+1;
				parent[nextunusednode-1] = tnode;
				parent[nextunusednode-1+1] = tnode;
				nextunusednode = nextunusednode+2;
			}
		} 

		if (leftchildren[tnode-1] == 0){ 
            /* Leaf => store results falling here (don't store them everywhere to avoid O(N^2) storage)*/
/*            printf("Leaf\n");*/

            /*=== Set up the structures to fill the cell array output. */
            dims[0] = Nnode;
            dims[1] = 1;
					
            mx_ysub = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
            mxSetCell(ysub, tnode-1, mx_ysub);
            ysub_for_result = mxGetPr(mx_ysub);

            mx_censsub = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
            mxSetCell(censsub, tnode-1, mx_censsub);
            censsub_for_result = (int*) mxGetData(mx_censsub);
            
            if (log10transformed){
/*                printf("log10\n");*/
                for(i=0; i<Nnode; i++){
                    ysub_for_result[i] = pow(10, ynode[i]);
                }
            } else {
/*                printf("no log10\n");*/
                for(i=0; i<Nnode; i++){
                    ysub_for_result[i] = ynode[i];
                }
            }
            for(i=0; i<Nnode; i++){
/*                if (log10transformed){ 
                    ysub_for_result[i] = pow(10, ynode[i]);
                } else { */
                ysub_for_result[i] = ynode[i];
/*                }*/
                censsub_for_result[i] = censnode[i];
            }
            
/*                printf("ynode[i]=%lf, ysub_for_result[%d]=%lf\n", ynode[i], i, ysub_for_result[i]);*/
        }
        
        tnode = tnode+1;
	}

	numNodesPointer[0] = nextunusednode - 1;
	numNcatsplitPointer[0] = ncatsplit;
	
	/*=== Free Memory for Auxiliary Arrays. */
    mxFree(noderows);
    mxFree(leftside);
    mxFree(rightside);
    mxFree(assignednode);
    mxFree(idx);
    mxFree(rows);
    mxFree(bestleft);
    mxFree(bestright);
    mxFree(censnode);
    mxFree(censnode_idx);
    mxFree(period);
    mxFree(period_idx);
    mxFree(N_all);
    mxFree(C_all);
    mxFree(O_all);

    mxFree(randomPermutation);
    
	mxFree(xnoderow);
	mxFree(ynode);
	mxFree(ynode_idx);
	mxFree(x);
    mxFree(period_time);
    
	mxFree(xleft);
	mxFree(xright);
	mxFree(numLeftPointer);
	mxFree(numRightPointer);
	mxFree(critvalPointer);
	mxFree(cutvalPointer);
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	double *X, *y, *tmp_double_ptr, *cutpoint;
	double p, percentageFeatures, kappa, cutoff_penalty_factor;
	int *tmp_int_ptr, *iscat, *nodenumber, *cens, *parent, *cutvar, *leftchildren, *rightchildren, *nodesize, *numNodesPointer, *numNcatsplitPointer, dim[1], dims[2];
	int numFeaturesType1, N, nvars, mrows, ncols, seed, SplitMinUncens, log10transformed, i;
	

  /* Check for proper number of arguments. */
  if(nrhs!=12 || nlhs != 12) {
    mexErrMsgTxt("USAGE: [nodenumber, parent, ysub, censsub, cutvar, cutpoint, leftchildren, rightchildren, nodesize, catsplit, numNodes, ncatsplit] = fh_random_regtreefit_logrank_big_leaves_twofeaturetypes(X, y, cens, SplitMinUncens, numFeaturesType1, p, percentageFeatures, iscat, kappa, cutoff_penalty_factor, seed, log10transformed).");
  }
  
  /* Check each argument for proper form and dimensions. */
  N = mxGetM(prhs[0]);
  nvars = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("X must be a noncomplex double matrix.");
  }
  X = mxGetPr(prhs[0]);

  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || !(mrows==N) || !(ncols==1) ) {
      printf("mrows=%d, N=%d\n",mrows,N);
    mexErrMsgTxt("y must be a noncomplex double column vector of the same length as size(X,1).");
  }
  y = mxGetPr(prhs[1]);

  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsInt32(prhs[2]) || mxIsComplex(prhs[2]) || !(mrows==N) || !(ncols==1) ) {
    mexErrMsgTxt("cens must be a noncomplex int column vector of the same length as size(X,1).");
  }
  cens = (int*) mxGetData(prhs[2]);
  
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsInt32(prhs[3]) || mxIsComplex(prhs[3]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("SplitMinUncens must be a noncomplex int scalar (cast it to int!).");
  }
  tmp_int_ptr = (int*) mxGetPr(prhs[3]);
  SplitMinUncens = tmp_int_ptr[0];
  
  
  mrows = mxGetM(prhs[4]);
  ncols = mxGetN(prhs[4]);
  if( !mxIsInt32(prhs[4]) || mxIsComplex(prhs[4]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("numFeaturesType1 must be a noncomplex int scalar (cast it to int!).");
  }
  tmp_int_ptr = (int*) mxGetPr(prhs[4]);
  numFeaturesType1 = tmp_int_ptr[0];

  mrows = mxGetM(prhs[5]);
  ncols = mxGetN(prhs[5]);
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("p must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[5]);
  p = tmp_double_ptr[0];

  mrows = mxGetM(prhs[6]);
  ncols = mxGetN(prhs[6]);
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("percentageFeatures must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[6]);
  percentageFeatures = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[7]);
  ncols = mxGetN(prhs[7]);
  if( !mxIsInt32(prhs[7]) || mxIsComplex(prhs[7]) || !(mrows==nvars) || !(ncols==1) ) {
    mexErrMsgTxt("iscat must be a noncomplex int column vector of the same length as size(X,2).");
  }
  iscat = (int*) mxGetPr(prhs[7]);

  mrows = mxGetM(prhs[8]);
  ncols = mxGetN(prhs[8]);
  if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("kappa must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[8]);
  kappa = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[9]);
  ncols = mxGetN(prhs[9]);
  if( !mxIsDouble(prhs[9]) || mxIsComplex(prhs[9]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("cutoff_penalty_factor must be a noncomplex double scalar.");
  }
  tmp_double_ptr = mxGetPr(prhs[9]);
  cutoff_penalty_factor = tmp_double_ptr[0];
  
  mrows = mxGetM(prhs[10]);
  ncols = mxGetN(prhs[10]);
  if( !mxIsInt32(prhs[10]) || mxIsComplex(prhs[10]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("seed must be a noncomplex int scalar.");
  }
  tmp_int_ptr = (int*) mxGetData(prhs[10]);
  seed = tmp_int_ptr[0];
  srand ( seed );
  
  mrows = mxGetM(prhs[11]);
  ncols = mxGetN(prhs[11]);
  if( !mxIsInt32(prhs[11]) || mxIsComplex(prhs[11]) || !(mrows==1) || !(ncols==1) ) {
    mexErrMsgTxt("log10transformed must be a noncomplex int scalar.");
  }
  tmp_int_ptr = (int*) mxGetData(prhs[11]);
  log10transformed = tmp_int_ptr[0];
  
  
  /* Create vectors for return arguments and assign pointers. */
  /* These have to be of size 2*N since the number of nodes can be that big (well, 2N-1) */
  dim[0] = 2*N;
  plhs[0] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  nodenumber = (int*) mxGetData(plhs[0]); 

  plhs[1] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  parent = (int*) mxGetPr(plhs[1]);

  plhs[2] = mxCreateCellArray(1, dim); /*ysub*/

  plhs[3] = mxCreateCellArray(1, dim); /*censsub*/

  plhs[4] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  cutvar = (int*) mxGetData(plhs[4]);

  plhs[5] = mxCreateNumericArray(1, dim, mxDOUBLE_CLASS, mxREAL);
  cutpoint = mxGetPr(plhs[5]);

  plhs[6] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  leftchildren = (int*) mxGetData(plhs[6]);

  plhs[7] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  rightchildren = (int*) mxGetData(plhs[7]);
  
  plhs[8] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  nodesize = (int*) mxGetData(plhs[8]);

  dims[0] = N;
  dims[1] = 2;
  plhs[9] = mxCreateCellArray(2, dims);

  dim[0] = 1;
  plhs[10] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  numNodesPointer = (int*) mxGetData(plhs[10]);

  plhs[11] = mxCreateNumericArray(1, dim, mxINT32_CLASS, mxREAL);
  numNcatsplitPointer = (int*) mxGetData(plhs[11]);
  
  /*=== Objective meanX, where X=cutoff_penalty_factor */
  /* We use this information directly to *build* the tree, so the tree
   * is already built using the metric we care about */
/*  printf("kappa = %lf\n", kappa);
  printMatrixDouble(y, 1, N);*/
  
/*
 if(log10transformed){
      kappa = log10(kappa);
  }
*/  
  
/*  printMatrixDouble(adjusted_y, 1, N); */
  
  /* printf("calling build the tree ...\n");  */
  
  buildTheTree(X, y, cens, SplitMinUncens, numFeaturesType1, p, percentageFeatures, iscat, N, nvars, log10transformed, kappa, cutoff_penalty_factor, nodenumber, parent, plhs[2], plhs[3], cutvar, cutpoint, leftchildren, rightchildren, nodesize, plhs[9], numNodesPointer, numNcatsplitPointer);
}