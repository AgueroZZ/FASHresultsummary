#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); //response variable
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  
  int d_beta = X.cols(); // Number of fixed effect coefficients

  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  DATA_SCALAR(intercept); // The assumed intercept
  // Parameter
  PARAMETER_VECTOR(W);
  int Wdim = W.size();
  int betadim = d_beta;
  vector<Type> beta(betadim);
  for (int i=0;i<betadim;i++) beta(i) = W(i);

  // Transformations
  vector<Type> eta = X * beta + intercept;

  // Log likelihood
  Type ll = 0;
  ll = sum(dpois(y, exp(eta), TRUE));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  lpW += -0.5 * d_beta * log(2*M_PI);

  REPORT(lpW);
  
  // Final result!
  Type logpost = -1 * (ll + lpW);
  
  return logpost;
}