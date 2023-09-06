
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
mat bayesian_lasso_rw_mode_c(mat mX, vec vy, double a, double b, double u, double v, 
                     mat mBeta, vec vsigma2, vec vlambda2)
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  const int maxiter = mBeta.n_rows;
  
  mat mMode(maxiter,p);
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  
  vec dgXTX  = (mX%mX).t() * one_n;
  vec XTy = mX.t() * vy;
  
  for (int i = 0; i < maxiter; ++i) 
  {
    ////////////////////////////////////////////////////////////////////////////

    vec va_val = dgXTX/vsigma2[i];
    double c_val = sqrt(vlambda2[i]/vsigma2[i]);
    
    vec vbeta = mBeta.row(i).as_col();
    vec vy_hat = mX * vbeta;
    
    for (int j=0; j<p; ++j)
    {
      vec vx_j = mX.col(j);
      vec vy_hat_mj = vy_hat - vx_j*vbeta[j];
      double b_val = (XTy[j] - sum(vx_j%vy_hat_mj))/vsigma2[i];
      mMode(i,j) = mlasso_c( va_val[j], b_val, c_val);
      vy_hat = vy_hat_mj + vx_j*vbeta[j];
    }
  }
  
  return mMode;
}




