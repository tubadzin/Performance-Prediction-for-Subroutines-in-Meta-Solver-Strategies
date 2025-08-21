/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>
#include <ctype.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include "stopwatch.h"
#include <set>

#include <vector>
#include <sstream>
#include <strings.h>
#include <stdio.h>
#include <algorithm>

using namespace std;
double hugeValue=1e100;
vector<string> featname;
vector<double> featvalue;
vector<string> timefeatname;
vector<double> timefeatvalue;
Stopwatch sw, sw2;
int num_stats_of_vector = 2;
double *tmp_stats = new double[num_stats_of_vector];
string *tmp_stats_names = new string[num_stats_of_vector];
int cur_numcols;
int cur_numrows;
int nzcnt;
double *obj_fun_coef;
double *rhs_b;
int *rmatbeg; 
int *rmatind; 
double *rmatval;
double noone=1234.1234;
string standardFeatString = "avg-median-varcoef-q90mq10";
double defaultForMissingValue = -512.0;

void strip(string &s);
void display(vector<double> &sv);
void getinfo(string mystring, string targetstr, vector<string> &sv);
void getinfosubstr(string mystring, int ss, int ee,  vector<string> &sv);
void str2double(vector<string> &sv, vector<double> &dv, int refill, double firstline, double bigvalue);
void splitVect(vector<string> &sv, vector<string> &sv1, vector<string> &sv2);
void splitoneVect(vector<string> &sv, vector<string> &sv1, vector<string> &sv2);
double str_to_f(string newstr);
void writefeat(string name, double value);
void writetimefeat(string name, double value);
void printfeat();
void printtimefeat();
void processfile(char* filename);
void basicVect(vector<double> &dv, double notcount, string feat_name, string which_statistics);
void countitemVect(vector<double>integralDV, double item, double &counter);     
void prob_cplex(CPXENVptr env, CPXLPptr lp);
double readfirst(vector<string> dv, double valueForEmpty);
void writeMfeature(string name, string flist, double mysum, double mymean, double mymax, double mymin, double mystd, double myq25, double mymedian, double myq75, double myvarcoef, double myinvvarcoef);
void diffVect(vector<double> dv, vector<double> &dv2, double notcount);
void computeFeatures_O_of_N(vector<double> array, double notcount, double &mysum, double &mymean,  double &mystd, double &myvarcoef, double &myinvvarcoef);
void computeFeatures_O_of_NlogN(vector<double> array, double notcount, double &mymax, double &mymin, double &myq10, double &myq25, double &mymedian, double &myq75, double &myq90 );
