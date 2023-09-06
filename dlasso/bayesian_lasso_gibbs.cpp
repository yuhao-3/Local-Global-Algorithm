
 
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <rinvgaussian.h>

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
vec rrinvgauss(vec vmu, vec vlambda)
{
  int n = vmu.n_rows;
  vec random_vector(n);
  double z,y,x,u;
  
  //Function f("rinvgauss"); 
  
  for (int i=0; i<n; ++i)
  {
    //random_vector[i] = f(1, _("mean")=vmu[i], _["disp"]=vlambda[i]));
    
    z = R::rnorm(0,1);
    y = z*z;
    x = vmu[i] + 0.5*vmu[i]*vmu[i]*y/vlambda[i] - 0.5*(vmu[i]/vlambda[i])*sqrt(4*vmu[i]*vlambda[i]*y + vmu[i]*vmu[i]*y*y);
    u = R::runif(0,1);
    if (u <= vmu[i]/(vmu[i]+x)) {
      random_vector[i] = x;
    } else {
      random_vector[i] = vmu[i]*vmu[i]/x;
    };
  }
  return(random_vector);
}

////////////////////////////////////////////////////////////////////////////////

// bayesian_lasso_gibbs_c
// Gibbs sampling for Bayesian lasso regression
// Assumes that lambda is fixed
// Assumes that n>p

// [[Rcpp::export]]
List bayesian_lasso_gibbs_c(mat mX, vec vy, double lambda, double sigma2_hat, 
                            double a, double b, int nburn, int nsamples, 
                            bool verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  mat XTX = mX.t()*mX;
  vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  int maxiter = nburn + nsamples;
  
  // Initialise storage of samples
  mat mBeta(maxiter,p);
  vec vsigma2(maxiter);
  
  // Initialize storage for Rao-Blackwellization
  mat mM(maxiter,p);
  mat mV(maxiter,p);
  vec va_til(maxiter);
  vec vb_til(maxiter);
  
  // Initialisation
  vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_hat;
  double lambda2 = lambda*lambda;
  double a_til = a + 0.5*(n + p);
  double b_til;
  
  mat vmu_til;
  mat mSigma_til;
  vec vbeta;
  mat mQ;
  
  // Parameters of va|rest which is inverse Gaussian
  vec vmu;
  vec vlambda = ones(p);
  
  
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    mQ = inv(XTX + diagmat(lambda2*va));
    vmu_til = mQ*XTy;
    mSigma_til = sigma2*mQ;
    vbeta = mvnrnd(vmu_til, mSigma_til, 1);
    
    // Sample from sigma2|rest
    b_til = b + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
   
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rrinvgauss(vmu, vlambda);
 
 
    if (verbose) {
      if ((i%1000)==1) {
        Rcout << "iter: " << i << "\n";
      }
    }
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = diagvec(mSigma_til).as_row();
    va_til[i] = a_til;
    vb_til[i] = b_til;
  }
  
  return List::create(_["mBeta"] = mBeta, 
                      _["vsigma2"] = vsigma2, 
                      _["mM"] = mM, 
                      _["mV"] = mV,
                      _["va_til"] = va_til, 
                      _["vb_til"] = vb_til);
}

////////////////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List bayesian_lasso_gibbs_tune_c(mat mX, vec vy, double lambda, double sigma2_hat, 
                            double a, double b, double u, double v,
                            int nburn, int nsamples, bool verbose) {
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  mat XTX = mX.t()*mX;
  vec XTy = mX.t()*vy;
  //double yTy = vy.t()*vy;
  
  int maxiter = nburn + nsamples;
  
  // Initialise storage of samples
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
  
  // Initialisation
  vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_hat;
  double lambda2 = lambda*lambda;
  double a_til = a + 0.5*(n + p);
  double b_til;
  double u_til = u + 0.5*p;
  double v_til;
  
  mat vmu_til;
  mat mSigma_til;
  vec vbeta;
  mat mQ;
  
  // Parameters of va|rest which is inverse Gaussian
  vec vmu;
  vec vlambda = ones(p);
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    mQ = inv(XTX + diagmat(lambda2*va));
    vmu_til = mQ*XTy;
    mSigma_til = sigma2*mQ;
    vbeta = mvnrnd(vmu_til, mSigma_til, 1);
    
    // Sample from sigma2|rest
    b_til = b + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(va%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rrinvgauss(vmu, vlambda);
    
    if (verbose) {
      if ((i%1000)==1) {
        Rcout << "iter: " << i << "\n";
      }
    }
    
    // Storing samples
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    mM.row(i) = vmu_til.as_row();
    mV.row(i) = diagvec(mSigma_til).as_row();
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


// looping through each column and element wise multiplication
// [[Rcpp::export]]
arma::mat matTimesVec(arma::mat mat, arma::vec v) {
  for(int i; i < mat.n_cols; i++){
    mat.col(i)  %=  v;
  }
  return mat;
}

// form a diagonal matrix with the vector and then use matrix multiplication
// [[Rcpp::export]]
arma::mat matTimesVec2(arma::mat mat, arma::vec v) {
  return arma::diagmat(v) * mat;
}

// use the functionality described at http://arma.sourceforge.net/docs.html#each_colrow 
// to "Apply a vector operation to each column or row of a matrix "
// [[Rcpp::export]]
arma::mat matTimesVec3(arma::mat mat, arma::vec v) {
  mat.each_col() %= v;
  return mat; 
}

mat crossprod_cpp(mat mX, vec vw, mat mD, double trunc) 
{
  int n = mX.n_rows;
  int p = mX.n_cols;
  
  mat out = mD;
  for (int i = 0; i<p; ++i) {
    if (vw[i]>trunc) {
      vec vx = mX.col(i);
      out = out + (vw[i]*vx) * vx.t();
    }
  }
  return out;
}




// [[Rcpp::export]]
List bayesian_lasso_gibbs_tune_pgtn_c(const mat& mX, const vec& vy, double lambda, 
                                      double sigma2_hat, const double a, const double b, 
                                      const double u, const double v, const int nburn, 
                                      const int nsamples, const int verbose, const double trunc) {
  const int n = mX.n_rows;
  const int p = mX.n_cols;
  
  // Summary statistics that can be calculated once at the beginning
  const mat& mI_n = eye(n, n);
  
  const int maxiter = nburn + nsamples;
  
  // Initialise storage of samples
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
  
  // Initialisation
  vec va = ones(p); // expected value of auxiliary variables under q
  double sigma2 = sigma2_hat;
  double lambda2 = lambda*lambda;
  const double a_til = a + 0.5*(n + p);
  double b_til;
  const double u_til = u + 0.5*p;
  double v_til;
  
  mat vmu_til;
  mat vsigma2_til;
  vec vbeta;
  mat mQ;
  mat mQ_inv;
  mat mC;
  mat XDh;
  //mat XYV = join_horiz(join_horiz(mX,vy),vy);
  mat QiXYV;
  mat XTQ;
  mat QX;
  mat XD;
  mat QXD;
  mat mU;
  mat mG;
  
  vec vd_inv;
  vec vd_invsqrt;
  vec vu;
  vec vv;
  vec vw;
  
  //Rcpp::Clock clock;
 
  // Parameters of va|rest which is inverse Gaussian
  vec vmu;
  const vec vlambda = ones(p);
  
  // Main loop
  for (int i = 0; i < maxiter; ++i) 
  {
    // Sample from beta|rest
    
    /**
    //clock.tick("A");
    vd_inv = 1/(lambda2*va);
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = -(mX*vu + randn(n));

    XD = mX.each_row() % vd_inv.t();
    
    if (trunc>0) {
      uvec ids = find( vd_inv>trunc );
      mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
    } else {
      mQ_inv = XD*mX.t() + mI_n;
    }
    
    mQ = inv_sympd(mQ_inv, inv_opts::allow_approx );
    QXD = mQ*XD;
    vmu_til = QXD.t()*vy;
    vw = (sum(XD % QXD,0)).t();
    vsigma2_til = sigma2*(vd_inv - vw);
    vw = QXD.t()*vv;
    vbeta = vmu_til + sqrt(sigma2)*(vu + vw);
    **/
    
    vd_inv = 1/(lambda2*va);
    vd_invsqrt = sqrt(vd_inv);
    vu = randn(p) % vd_invsqrt;
    vv = mX*vu + randn(n);
    
    XD = mX.each_row() % vd_inv.t();
    if (trunc>0) {
      uvec ids = find( vd_inv>trunc );
      mQ_inv = XD.cols(ids)*mX.cols(ids).t() + mI_n;
    } else {
      mQ_inv = XD*mX.t() + mI_n;
    }
    
    // Force symmetry
    // Sometimes for numerical reasons it isn't symmetric.
    mQ_inv = 0.5*(mQ_inv + mQ_inv.t());
    
    mU = trimatu(chol(mQ_inv));
    mat UT = trimatl(mU.t()); 
        
    mG = solve(UT, mX);
    vw = (mG % mG).t() * one_n; 
    vsigma2_til = sigma2*(vd_inv - vd_inv%vd_inv%vw);
        
    vw = solve(mU, solve(mU.t(), -vv));
    vmu_til = mX.t() * solve(mU, solve(UT, vy));
    vmu_til = vd_inv % vmu_til;
    vbeta = vmu_til + sqrt(sigma2)*(vu + vd_inv%(mX.t()*vw));

    
    // Sample from sigma2|rest
    b_til = b + 0.5*(sum(pow(vy-mX*vbeta,2.0)) + lambda2*sum(va%pow(vbeta,2.0)));
    sigma2 = 1/randg(distr_param(a_til,1/b_til));
    
    // Sample from lambda2|rest
    v_til = v + 0.5*sum(va%pow(vbeta,2.0))/sigma2;
    lambda2 = randg(distr_param(u_til,1/v_til));
    
    // Update q(va) 
    vmu = sqrt(sigma2/(lambda2*(pow(vbeta,2.0))));
    va = rrinvgauss(vmu, vlambda);
    
    if (verbose!=0) {
      if ((i%verbose)==0) {
        Rcout << "iter: " << i << " lambda2: " << lambda2 << " sigma2: " << sigma2 << "\n";
      }
    }
    
    // Storing samples
    //Rcout << "A: \n";
    mBeta.row(i) = vbeta.as_row();
    vsigma2[i] = sigma2;
    vlambda2[i] = lambda2;
    
    // Storage for Rao-Blackwellization
    //Rcout << "B: \n";
    mM.row(i) = vmu_til.as_row();
    //Rcout << "C: \n";
    mV.row(i) = vsigma2_til.as_row();
    //Rcout << "D: \n";
    va_til[i] = a_til;
    vb_til[i] = b_til;
    vu_til[i] = u_til;
    vv_til[i] = v_til;
  }
  
  //clock.stop("Rcpp_times");
  
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

