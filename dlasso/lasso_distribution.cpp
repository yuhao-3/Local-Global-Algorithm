

#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "zeta.h"
 
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

typedef Map<MatrixXd> MapMat;
typedef Map<VectorXd> MapVec;

////////////////////////////////////////////////////////////////////////////////

#include <R.h>
#include <Rmath.h>
 
/* Global Variables */
#define TINY1   1.0E-30
#define TINY2   1.0E-7

////////////////////////////////////////////////////////////////////////////////

// The inverse logit function

// [[Rcpp::export]]
double expit_c(double x) {
  double val = 1/(1 + exp(-x));
  return val;
}

////////////////////////////////////////////////////////////////////////////////

Rcpp::List calculate_lasso_dist_stats_c(double a_val, double b_val, double c_val) 
{
  double mu_plus = (b_val - c_val)/a_val;
  double mu_minus = (b_val + c_val)/a_val;
  double sigma2 = 1/a_val;
  double sigma = sqrt(sigma2);
  
  double r_plus = mu_plus/sigma;
  double r_minus = mu_minus/sigma;
  
  double z_plus  = zeta_c(1, r_plus);
  double z_minus = zeta_c(1,-r_minus);
  
  double log_pm = R::pnorm(-(b_val + c_val)*sigma, 0.0, 1.0, 1, 1);
  double log_pp = R::pnorm( (b_val - c_val)*sigma, 0.0, 1.0, 1, 1);
  double w = expit_c(log_pm - log_pp + 2*b_val*c_val*sigma2);
  
  return List::create(_["mu_plus"] = mu_plus, 
                      _["mu_minus"] = mu_minus, 
                      _["sigma2"] = sigma2, 
                      _["sigma"] = sigma, 
                      _["z_plus"] = z_plus, 
                      _["z_minus"] = z_minus,
                      _["r_plus"] = r_plus, 
                      _["r_minus"] = r_minus, 
                      _["w"] = w);
}

////////////////////////////////////////////////////////////////////////////////

double logSumExp_c(vec vx) {
  double M = max(vx);
  double val = M + log(sum(exp(vx - M)));
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// return the normalizing constant

// [[Rcpp::export]]
double zlasso_c(double a_val, double b_val, double c_val, bool logarithm) 
{
  List res = calculate_lasso_dist_stats_c(a_val, b_val, c_val);
  double r_plus = res["r_plus"];
  double r_minus = res["r_minus"];
  double sigma = res["sigma"];
  
  double log_inv_z_plus  = R::pnorm5( r_plus, 0.0, 1.0, 1, 1) - R::dnorm4( r_plus, 0.0, 1.0, 1);
  double log_inv_z_minus = R::pnorm5(-r_minus,0.0, 1.0, 1, 1) - R::dnorm4(-r_minus,0.0, 1.0, 1);
  
  vec vx(2);
  vx[0] = log_inv_z_plus;
  vx[1] = log_inv_z_minus;
  double val = log(sigma) + logSumExp_c(vx);
  if (logarithm) {
    return val;
  } 
  val = exp(val);
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// Note: a>0, c>0

// [[Rcpp::export]]
double dlasso_c(double x, double a_val, double b_val, double c_val, bool logarithm) 
{
  double log_Z = zlasso_c(a_val, b_val, c_val, true); 
  double val =  -0.5*a_val*x*x + b_val*x - c_val*abs(x) - log_Z;
  
  if (logarithm) {
    return val;  
  }  
  val = exp(val);
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double plasso_c(double x, double a_val, double b_val, double c_val) 
{
  List res = calculate_lasso_dist_stats_c(a_val, b_val, c_val);
  double w = res["w"];
  double mu_plus = res["mu_plus"];
  double mu_minus = res["mu_minus"];
  double r_minus = res["r_minus"];
  double r_plus = res["r_plus"];
  double sigma = res["sigma"];
  
  double val = 0.0;
  if (x<=0) {
    val = w*exp( R::pnorm5((x-mu_minus)/sigma, 0.0, 1.0, 1, 1) - R::pnorm5(-r_minus, 0.0, 1.0, 1, 1));
  } else {
    val = w + (1.0 - w)*(1.0 - exp( R::pnorm5((mu_plus - x)/sigma, 0.0, 1.0, 1, 1) - R::pnorm5(r_plus, 0.0, 1.0, 1, 1)));
  }
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double qlasso_fast_c(double u, double a_val, double b_val, double c_val) 
{
  double sigma2 = 1/a_val;
  double sigma = sqrt(sigma2);
  
  double log_pm = R::pnorm5(-(b_val + c_val)*sigma, 0.0, 1.0, 1, 1);
  double log_pp = R::pnorm5( (b_val - c_val)*sigma, 0.0, 1.0, 1, 1);
  double w = expit_c(log_pm - log_pp + 2*b_val*c_val*sigma2);
  
  double x;
  if (u<=w) {
    x =  (b_val + c_val)*sigma2 + sigma*R::qnorm5( exp(log_pm)*u/w, 0.0, 1.0, 1, 0);
  } else {
    x =  (b_val - c_val)*sigma2 - sigma*R::qnorm5( exp(log_pp)*(1 - u)/(1-w), 0.0, 1.0, 1, 0);
  }
  return x;
}

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double rlasso_fast_c(double a_val, double b_val, double c_val) {
  double u = R::runif(0.0, 1.0);
  double x = qlasso_fast_c(u, a_val, b_val, c_val);
  return x;
}

////////////////////////////////////////////////////////////////////////////////

// return the expected value

// [[Rcpp::export]]
double elasso_c(double a_val, double b_val, double c_val) 
{
  List res = calculate_lasso_dist_stats_c(a_val, b_val, c_val);
  double r_plus = res["r_plus"];
  double r_minus = res["r_minus"];
  double mu_plus = res["mu_plus"];
  double mu_minus = res["mu_minus"];
  double sigma = res["sigma"];
  double w = res["w"];
  
  double z_plus = zeta_c(1, r_plus);
  double z_minus = zeta_c(1,-r_minus);
  double e_plus  =  mu_plus  + sigma*z_plus;
  double e_minus =  mu_minus - sigma*z_minus;
  double val = w*e_minus + (1.0-w)*e_plus; 
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// return the variance

// [[Rcpp::export]]
double vlasso_c(double a_val, double b_val, double c_val) 
{
  List res = calculate_lasso_dist_stats_c(a_val, b_val, c_val);
  double r_plus = res["r_plus"];
  double r_minus = res["r_minus"];
  double mu_plus = res["mu_plus"];
  double mu_minus = res["mu_minus"];
  double sigma = res["sigma"];
  double sigma2 = res["sigma2"];
  double w = res["w"];
  
  double z_plus = zeta_c(1, r_plus);
  double z_minus = zeta_c(1,-r_minus);
  double e_plus = mu_plus  + sigma*z_plus;
  double e_minus = mu_minus - sigma*z_minus;
  
  double v_plus  = sigma2*(1 + zeta_c(2, r_plus));
  double v_minus = sigma2*(1 + zeta_c(2,-r_minus));
  
  double val = w*(v_minus + e_minus*e_minus)  + (1.0 - w)*(v_plus + e_plus*e_plus);
  val = val - pow( w*e_minus + (1-w)*e_plus, 2.0);
  
  return val;
}

////////////////////////////////////////////////////////////////////////////////

// Return the mode of the lasso distribution

// [[Rcpp::export]]
double mlasso_c(double a_val, double b_val, double c_val) {
  vec vx(2);
  vx[0] = abs(b_val) - c_val;
  vx[1] = 0.0;
  double val = max(vx)*sign(b_val)/a_val;
  return val;
}

////////////////////////////////////////////////////////////////////////////////
