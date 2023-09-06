
#ifndef LASSO_DISTRIBUTION_H
#define LASSO_DISTRIBUTION_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////////////////////////////////////

double expit_c(double x);
//Rcpp::List calculate_lasso_dist_stats_c(double a_val, double b_val, double c_val);
//double logSumExp_c(vec vx);
double zlasso_c(double a_val, double b_val, double c_val, bool logarithm);
double dlasso_c(double x, double a_val, double b_val, double c_val, bool logarithm);
double plasso_c(double x, double a_val, double b_val, double c_val);
double qlasso_fast_c(double u, double a_val, double b_val, double c_val);
double rlasso_fast_c(double a_val, double b_val, double c_val);
double elasso_c(double a_val, double b_val, double c_val);
double vlasso_c(double a_val, double b_val, double c_val);
double mlasso_c(double a_val, double b_val, double c_val);

////////////////////////////////////////////////////////////////////////////////

#endif
