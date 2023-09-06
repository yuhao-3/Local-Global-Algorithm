
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"
#include "slice_sampler.h"
#include "lasso_distribution.h"

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


// [[Rcpp::export]]
List lasso_gibbs_Bc(mat mX, vec vy, double a, double b, double u, double v, 
                     int nsamples, vec vbeta_init, double lambda_init, double sigma2_init,
                     int verbose)
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples for MCMC
  mat mBeta(maxiter,p);
  vec vsigma2(maxiter);
  vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  mat mA(maxiter,p);
  mat mB(maxiter,p);
  mat mC(maxiter,p);
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  mat XTX;
  vec dgXTX  = (mX%mX).t() * one_n;
  if (n>p) {
    XTX = mX.t() * mX;
  } 
  vec XTy = mX.t() * vy;
  
  //Rcout << "A: \n";
  
  // Set the current values of the parameters 
  vec vb = ones(p); // expected value of auxiliary variables under q
  vec vnu = ones(p);
  
  // Assign initial values
  vec vbeta = vbeta_init;
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda_init*lambda_init;
  
  // Constant values 
  const double a_til = a + 0.5*(n + p);
  double b_til = b;
  const double u_til = u + 0.5*p;
  double v_til = v;
  
  //Rcout << "B: \n";
  
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from vb|rest
    vnu = sigma/(lambda*abs(vbeta));  
    vb  = rinvgaussian_c(vnu, one_p);  
    
    double num;
    double denom;
    
    vec va_vals = dgXTX/sigma2;
    vec vb_vals = zeros(p);
    vec vc_vals = one_p*sqrt(lambda2/sigma2);
    
    vec vu = randu(p);
    if (n>p) {
      vec XTy_hat = XTX * vbeta;
      for (int j=0; j<p; ++j) {
        vec vx_j = XTX.col(j);
        XTy_hat = XTy_hat - vx_j*vbeta[j]; // This might not be exactly right
        num = XTy[j] - XTy_hat[j];
        vb_vals[j] = num/sigma2;
        denom = dgXTX[j] + lambda2*vb[j];
        vbeta[j] =  qlasso_fast_c(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        XTy_hat = XTy_hat + vx_j*vbeta[j];
      }
    } else {
      vec vy_hat = mX * vbeta;
      for (int j=0; j<p; ++j) {
        vec vx_j = mX.col(j);
        vec vy_hat_mj = vy_hat - vx_j*vbeta[j];
        num = XTy[j] -  as_scalar(vx_j.t() * vy_hat_mj);
        vb_vals[j] = num/sigma2;
        vbeta[j] = qlasso_fast_c(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
        vy_hat = vy_hat_mj +  vx_j*vbeta[j];
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Slice from lambda2|rest
    double sum_abs_vbeta = sum(abs(vbeta));
    double RSS  = sum(pow(vy-mX*vbeta,2.0));
    
    // Slice sample from sigma2|rest
    lambda = sqrt(lambda2);
      
    double a_val = (a_til-1);
    double b_val = b + 0.5*RSS;
    double c_val = lambda*sum_abs_vbeta;
    
    // Slice sampler for tau and then invert.
    double tau = 1/sigma2;
    tau = slice_sampler_precision_c(tau, a_val, b_val, c_val);
    sigma2 = 1/tau;
    sigma = sqrt(sigma2);
    
    ////////////////////////////////////////////////////////////////////////////
    
    a_val = u_til - 1;
    b_val = v;
    c_val = sum_abs_vbeta/sigma;
    
    // Slice sampler for lambda2
    lambda2 = slice_sampler_precision_c(lambda2, a_val, b_val, c_val);
    lambda  = sqrt(lambda2);
    
    ////////////////////////////////////////////////////////////////////////////
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    // Store MCMC samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage sufficient statistics for Rao-Blackwellization
    mA.row(i) = va_vals.as_row();
    mB.row(i) = vb_vals.as_row();
    mC.row(i) = vc_vals.as_row();
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mA"] = mA, 
                      _["mB"] = mB,
                      _["mC"] = mC);
}




