
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include "rinvgaussian.h"

//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>

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

// [[Rcpp::export]]
arma::vec subsum(const arma::mat& X, const arma::vec& vbeta, arma::uword i) {
  arma::uword n = X.n_rows; // X should be square so only need # rows
  arma::uvec idx(n-1); // vector of indices to subset by
  arma::uword ii = 0; // the integer we'll add at each elem of idx
  for ( arma::uword j = 0; j < (n-1); ++j ) { // for each elem of idx
    if ( ii == i ) { // if ii equals i, we need to skip i
      ii += 1;     // (i.e., add 1 to ii)
    }
    idx[j] = ii;     // then we store ii for this elem
    ii += 1;         // and increment ii
  }
  uvec jj(2);
  jj[0] = i;
  jj[1] = i;
  
  uvec oneone(2);
  jj[0] = 0;
  jj[1] = 0;
  arma::vec res =  X.submat(jj, idx) * vbeta.submat(idx,oneone);
  return res; 
}

/** 

vec rinvgauss_c(vec vmu, vec vlambda) 
{
  int n = vmu.n_rows;
  vec vphi = 1/vlambda;
  vec vr = zeros(n);  
  uvec cond1 = find((vmu> 0) && (vphi> 0));
  uvec cond2 = find((vmu<=0) || (vphi<=0));
  
  // Take care of samples with bad inputs
  if (cond2.n_elem>0) {
    vr.elem(cond2).fill(datum::nan);
    n = cond1.n_elem;
  }
  
  // For samples with good arguments calculate vy and vx
  vphi.elem(cond1) = vphi.elem(cond1) % vmu.elem(cond1);
  vec vy = chi2rnd( 1.0, n );
  vec vx = 1.0 + vphi.elem(cond1)/2.0 % (vy - sqrt(4*vy/vphi.elem(cond1) + vy%vy));
  
  // Note: The line above should yield all vx>0, but it occasionally doesn't due to
  // numerical precision issues. The line below detects this and recomputes
  // the relevant elements of vx using a 2nd-order Taylor expansion of the
  // sqrt function, which is a good approximation whenever the problem occurs.
  
  uvec cond3 = find(vx<=0);
  if (cond3.n_elem>0) {
    vx.elem(cond3) = (1/(vy.elem(cond3) % vphi.elem(cond3)));
  }
  
  vec vu = randu(n);
  uvec cond4 = find(vu< (1.0/(1+vx)));
  uvec cond5 = find(vu>=(1.0/(1+vx)));
  
  vec temp = vr.elem(cond1);
  temp.elem(cond4) = vx.elem(cond4);
  temp.elem(cond5) = 1/vx.elem(cond5);
  vr.elem(cond1) = temp;
  
  return (vmu%vr);
}

// Draw from inverse-Gaussian distribution while avoiding potential numerical problems


// [[Rcpp::export]]
vec rinvgaussian_c(vec vmu, vec vlambda) {
  vec vm = vmu / sqrt(vmu % vlambda);
  vec vl = vlambda / sqrt(vmu % vlambda);
  vec result = sqrt(vmu % vlambda) % rinvgauss_c(vm, vl);
  return result;
}

**/

////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List lasso_gibbs_A_c(mat mX, vec vy, double a, double b, double u, double v, 
                     int nsamples, vec vbeta_init, double lambda_init, double sigma2_init)
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  const int maxiter = nsamples;
  
  // Initialise storage of samples for MCMC
  mat mBeta(maxiter,p);
  vec vsigma2(maxiter);
  vec vlambda2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  mat mM(maxiter,p);
  mat mV(maxiter,p);
  vec va_til(maxiter);
  vec vb_til(maxiter);
  vec vu_til(maxiter);
  vec vv_til(maxiter);
  
  const vec one_n = ones(n);
  const vec one_p = ones(p);
  mat XTX;
  vec dgXTX  = (mX%mX).t() * one_n;
  if (n>p) {
    XTX = mX.t() * mX;
  } 
  vec XTy = mX.t() * vy;
  
  // Set the current values of the parameters 
  vec vb = ones(p); // expected value of auxiliary variables under q
  vec vnu = ones(p);
  vec vsigma2_til=ones(p);
  vec vsigma_til=ones(p);
  vec vmu_til=zeros(p);
  vec vbeta = vbeta_init;
  double sigma2 = sigma2_init;
  double sigma = sqrt(sigma2);
  double lambda = lambda_init;
  double lambda2 = lambda_init*lambda_init;
  const double a_til = a + 0.5*(n + p);
  double b_til = b;
  const double u_til = u + 0.5*p;
  double v_til = v;
  
  for (int i = 0; i < maxiter; ++i) 
  {
    sigma  = sqrt(sigma2);
    lambda = sqrt(lambda2);
    
    vnu = sigma/(lambda*abs(vbeta));  
    vb  = rinvgaussian_c(vnu, one_p);  
    
    double num;
    double denom;
    vsigma2_til = sigma2/(dgXTX + lambda2*vb);
    vsigma_til = sqrt(vsigma2_til);
    
    if (n>p) {
       vec XTy_hat = XTX * vbeta;
       for (int j=0; j<p; ++j) {
         vec vx_j = XTX.col(j);
         XTy_hat = XTy_hat - vx_j*vbeta[j]; // This might not be exactly right
         num = XTy[j] - XTy_hat[j];
         denom = dgXTX[j] + lambda2*vb[j];
         vmu_til[j] = num/denom;
         vbeta[j] = randn(distr_param(vmu_til[j], vsigma_til[j]));
         XTy_hat = XTy_hat + vx_j*vbeta[j];
       }
    } else {
       vec vy_hat = mX * vbeta;
       for (int j=0; j<p; ++j) 
       {
          vec vx_j = mX.col(j);
          vec vy_hat_mj = vy_hat - vx_j*vbeta[j];
          num = XTy[j] -  as_scalar(vx_j.t() * vy_hat_mj);
          denom = dgXTX[j] + lambda2*vb[j];
          vmu_til[j] = num/denom;
          vbeta[j] =  randn(distr_param(vmu_til[j], vsigma_til[j]));
          vy_hat = vy_hat_mj +  vx_j*vbeta[j];
       }
    }
     
    // Sample from sigma2|rest
    b_til = b + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(vb%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(vb%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));

    // Store MCMC samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage sufficient statistics for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = vsigma2_til.as_row();
    va_til[i] = a_til;
    vb_til[i] = b_til;
    vu_til[i] = u_til;
    vv_til[i] = v_til;
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["vlambda2"] = vlambda2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["va_til"] = va_til, 
                      _["vb_til"] = vb_til,
                      _["vu_til"] = vu_til, 
                      _["vv_til"] = vv_til);
}
  
