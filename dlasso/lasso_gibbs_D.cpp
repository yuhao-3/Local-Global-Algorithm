
#define ARMA_WARN_LEVEL 0

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
List lasso_gibbs_Dc(mat mX, vec vy, double a, double b, double u, double v, 
                     int nsamples, vec vbeta_init, double lambda_init, double sigma2_init,
                     vec sI, int verbose)
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  const int maxiter = nsamples;
  
  // Index sets of elements in sI and sJ
  vec sJ = 1.0 - sI;
  arma::uvec sI_ids = find(sI==1);
  arma::uvec sJ_ids = find(sJ==1);
  
  int s = (int)(sum(sI));
  int r = p - s;
  
  // Initialise storage of samples for MCMC
  mat mBeta(maxiter,p);
  vec vsigma2(maxiter);
  vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  mat mM(maxiter,s);
  mat mV(maxiter,s);
  mat mA(maxiter,r);
  mat mB(maxiter,r);
  mat mC(maxiter,r);
  mat mMode(maxiter,r);
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  const vec one_r = ones(r);
  const vec one_s = ones(s);
  
  mat XTX;
  vec dgXTX  = (mX%mX).t() * one_n;
  
  if (n>p) {
    XTX = mX.t() * mX;
  } 
  vec XTy = mX.t() * vy;
  
  
  mat mX1 = mX.cols(sI_ids);
  mat mX0 = mX.cols(sJ_ids);
  
  vec dgXTX1 = dgXTX.elem(sI_ids);
  vec dgXTX0 = dgXTX.elem(sJ_ids);
  
  mat X1TX1 = mX1.t() * mX1;
  mat X1TX0 = mX1.t() * mX0;
  vec X1Ty = mX1.t() * vy;
  vec X0Ty = mX0.t() * vy;
 
  // Set the current values of the parameters 
  vec vb = ones(s);  
  vec vnu = ones(s);
  
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

  double num;
  double denom;
  
  int j;
  
  mat mQ_inv;
  mat mQ;
  mat mSigma_til;
  vec vsigma2_til;
  vec vmu_til;
  
  for (int i = 0; i < maxiter; ++i) 
  {
    vec vbeta_sI = vbeta.elem(sI_ids);
    vec vbeta_sJ = vbeta.elem(sJ_ids);
    
    //vbeta_sI.print();
    
    // Sample from vb|rest
    vnu = sigma/(lambda*abs(vbeta_sI));  
    vb  = rinvgaussian_c(vnu, one_s); 
 
    ////////////////////////////////////////////////////////////////////////////
    
    // Sample from vbeta_sI|rest
    mQ_inv = X1TX1 + diagmat(lambda2*vb);
    mQ = inv(mQ_inv);
    mSigma_til = sigma2*mQ;
    vsigma2_til = diagvec(mSigma_til);
    vec vmu1 = mQ * X1Ty;
    mat mP = mQ * X1TX0;
    vmu_til = vmu1 - mP*vbeta_sJ;
    vbeta_sI = mvnrnd(vmu_til, mSigma_til, 1);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Fill in elements of vbeta
    for (int jj=0; jj<s; ++jj) {
      j = sI_ids[jj];
      vbeta[j] = vbeta_sI[jj];
    }

    ////////////////////////////////////////////////////////////////////////////

    vec vu = randu(r);
    vec va_vals;
    vec vb_vals = zeros(r);
    vec vc_vals = one_r*sqrt(lambda2/sigma2);
    
 
    va_vals = (dgXTX0 - (X1TX0 % mP).t()*one_s)/sigma2;
    mat mX0_hat = mX0 - mX1 * mP;
    vec vy_hat = mX1*vmu1 + mX0_hat * vbeta_sJ;
    
    for (int j=0; j<r; ++j)
    {
      vec vx_j_hat = mX0_hat.col(j);
      vec vx_j = mX0.col(j);
      
      vec vy_hat_mj = vy_hat - vx_j_hat*vbeta_sJ[j];
      vb_vals[j] = (X0Ty[j] - sum(vx_j%vy_hat_mj))/sigma2;
      vbeta_sJ[j] = qlasso_fast_c(vu[j], va_vals[j], vb_vals[j], vc_vals[j]);
      mMode(i,j) = mlasso_c(va_vals[j], vb_vals[j], vc_vals[j]);
      vy_hat = vy_hat_mj + vx_j_hat*vbeta_sJ[j];
    }
    
    vy_hat = mX1 * vbeta_sI + mX0 * vbeta_sJ;
    
    
    /**
    if (n>p) {
      vec XTy_hat = XTX * vbeta;
      for (int jj=0; jj<r; ++jj) {
        j = sJ_ids[jj];
        vec vx_j = XTX.col(j);
        XTy_hat = XTy_hat - vx_j*vbeta[j];  
        num = XTy[j] - XTy_hat[j];
        vb_vals[jj] = num/sigma2;
        denom = dgXTX[j] + lambda2*vb[j];
        vbeta[j] =  qlasso_fast_c(vu[jj], va_vals[jj], vb_vals[jj], vc_vals[jj]);
        XTy_hat = XTy_hat + vx_j*vbeta[j];
      }
    } else {
      vec vy_hat = mX * vbeta;
      for (int jj=0; jj<r; ++jj) {
        j = sJ_ids[jj];
        vec vx_j = mX.col(j);
        vec vy_hat_mj = vy_hat - vx_j*vbeta[j];
        num = XTy[j] -  as_scalar(vx_j.t() * vy_hat_mj);
        vb_vals[jj] = num/sigma2;
        vbeta[j] = qlasso_fast_c(vu[jj], va_vals[jj], vb_vals[jj], vc_vals[jj]);
        vy_hat = vy_hat_mj +  vx_j*vbeta[j];
      }
    }
     
    **/
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Fill in elements of vbeta
    for (int jj=0; jj<r; ++jj) {
      j = sJ_ids[jj];
      vbeta[j] = vbeta_sJ[jj];
    }
    
    
    //vbeta_sI = vbeta.elem(sI_ids);
    //vbeta_sJ = vbeta.elem(sJ_ids);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Calculate some shared values used for sampling from sigma2 and lambda2
    double sum_abs_vbeta = sum(abs(vbeta_sJ));
    double RSS  = sum(pow(vy-vy_hat,2.0));
    double sum_vbeta_sq_on_va = sum(vb % pow(vbeta_sI,2));
    
    //Rcout << "sum_abs_vbeta: " << sum_abs_vbeta << "\n";
    //Rcout << "RSS: " << RSS << "\n";
    //Rcout << "sum_vbeta_sq_on_va: " << sum_vbeta_sq_on_va << "\n";
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Slice sample from sigma2|rest
    lambda = sqrt(lambda2);
      
    double a_val = (a_til-1);
    double b_val = b + 0.5*RSS + 0.5*sum_vbeta_sq_on_va*lambda2;
    double c_val = lambda*sum_abs_vbeta;
    
    //Rcout << "Slice sample from sigma2|rest: " << "\n";
    //Rcout << "a_val: " << a_val << "\n";
    //Rcout << "b_val: " << b_val << "\n";
    //Rcout << "c_val: " << c_val << "\n";
    //Rcout << "lambda2: " << lambda2 << "\n";
    //Rcout << "sigma2: " << sigma2 << "\n";
    
    // Slice sampler for tau and then invert.
    double tau = 1/sigma2;
    tau = slice_sampler_precision_c(tau, a_val, b_val, c_val);
    sigma2 = 1/tau;
    
    //Rcout << "output: " << sigma2 << "\n";
    
    sigma = sqrt(sigma2);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Slice from lambda2|rest
    a_val = u_til - 1;
    b_val = v + 0.5*sum_vbeta_sq_on_va/sigma2;
    c_val = sum_abs_vbeta/sigma;
    
    //Rcout << "Slice sample from lambda2|rest: " << "\n";
    //Rcout << "a_val: " << a_val << "\n";
    //Rcout << "b_val: " << b_val << "\n";
    //Rcout << "c_val: " << c_val << "\n";
    //Rcout << "lambda2: " << lambda2 << "\n";
    //Rcout << "sigma2: " << sigma2 << "\n";
    
    // Slice sampler for lambda2
    lambda2 = slice_sampler_precision_c(lambda2, a_val, b_val, c_val);
    lambda  = sqrt(lambda2);
    
    //Rcout << "output: " << lambda2 << "\n";
    
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
    
    // Storage sufficient statistics for Rao-Blackwellization for sI
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = vsigma2_til.as_row();
    
    // Storage sufficient statistics for Rao-Blackwellization for sJ
    mA.row(i) = va_vals.as_row();
    mB.row(i) = vb_vals.as_row();
    mC.row(i) = vc_vals.as_row();
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["mA"] = mA, 
                      _["mB"] = mB,
                      _["mC"] = mC,
                      _["mMode"] = mMode);
}




