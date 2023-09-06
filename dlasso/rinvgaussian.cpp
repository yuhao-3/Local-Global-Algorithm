

#include <RcppArmadillo.h>
#include <RcppNumerical.h>

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
  
  for (int i=0; i<n; ++i)
  {
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

  /** 
   * Note: The line above should yield all vx>0, but it occasionally doesn't due to
   * numerical precision issues. The line below detects this and recomputes
   * the relevant elements of vx using a 2nd-order Taylor expansion of the
   * sqrt function, which is a good approximation whenever the problem occurs.
   **/
  
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

/**
Draw from inverse-Gaussian distribution while avoiding potential numerical problems
**/

// [[Rcpp::export]]
vec rinvgaussian_c(vec vmu, vec vlambda) {
  vec vm = vmu / sqrt(vmu % vlambda);
  vec vl = vlambda / sqrt(vmu % vlambda);
  vec result = sqrt(vmu % vlambda) % rinvgauss_c(vm, vl);
  return result;
}
