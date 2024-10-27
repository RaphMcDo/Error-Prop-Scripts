// ------------------------------------------------------------------------
// modeling of weight~height+age+... : spatial slope
// ------------------------------------------------------------------------

#include <TMB.hpp>

#include "distmat.hpp"

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() (){
  
  using namespace density;
  
  
  DATA_VECTOR(weight);
  // DATA_UPDATE(weight);
  DATA_MATRIX(varmat);
  // DATA_UPDATE(varmat);
  DATA_VECTOR(depth);
  // DATA_UPDATE(depth);
  DATA_FACTOR(ind_loc);
  // DATA_UPDATE(ind_loc);
  DATA_MATRIX(locations); 
  // DATA_UPDATE(locations);
  
  //Spatial params
  PARAMETER(log_rho);
  PARAMETER(log_sig);
  PARAMETER(log_nu);
  PARAMETER(log_rho_b);
  PARAMETER(log_sig_b);
  PARAMETER(log_nu_b);
  
  Type rho = exp(log_rho);
  Type sig = exp(log_sig);
  Type nu = exp(log_nu);
  Type rho_b = exp(log_rho_b);
  Type sig_b = exp(log_sig_b);
  Type nu_b = exp(log_nu_b);
  
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(beta_s);
  PARAMETER_VECTOR(beta_b_s);
  PARAMETER_VECTOR(beta_depth);
  PARAMETER(log_phi);
  
  
  int n_loc = beta_s.rows(); 
  // int n_var = varmat.cols();
  
  
  // ---------------------------
  // Joint negative log-likelihood
  Type jnll = Type(0.0);
  
  // Spatial effect: Random intercept and Random slope
  
  matrix <Type> D(n_loc,n_loc);
  D = DistMat2d(locations);
  matrix<Type> St(n_loc,n_loc);
  St.setZero();
  matrix<Type> St_b(n_loc,n_loc);
  St_b.setZero();
  for(int i=0; i<n_loc; ++i){
    St(i,i) = sig*sig;
    St_b(i,i) = sig_b*sig_b;
  }
  for(int i=0; i<n_loc; ++i){
    for(int j=i+1; j<n_loc; ++j){
      St(i,j) = sig*sig*matern(D(i,j), rho, nu); //Matern Covariance Function
      St(j,i) = St(i,j);
      St_b(i,j) = sig_b*sig_b*matern(D(i,j), rho_b, nu_b); //Matern Covariance Function
      St_b(j,i) = St_b(i,j);
    }
  }
  
  jnll += density::MVNORM(St)(beta_s);
  jnll += density::MVNORM(St_b)(beta_b_s);
  
  // SIMULATE{
  //   vector<Type>temp_beta_s(n_loc);
  //   density::MVNORM(St).simulate(temp_beta_s);
  //   beta_s = temp_beta_s;
  //   REPORT(beta_s);
  // }
  
  // Observation likelihood
  int nobs = weight.size();
  vector<Type> mu(nobs); mu.setZero(); 
  vector<Type> mu_fixed(nobs); mu_fixed.setZero(); 
  vector<Type> mu_s(nobs); mu_s.setZero();
  vector<Type> mu_depth(nobs); mu_depth.setZero();  
  for (int i = 0; i < nobs; i++){
    for(int j = 0; j < varmat.cols(); j++){
      if (j==1) mu_fixed(i) += varmat(i, j) * (beta(j)+beta_b_s(ind_loc(i)));
      if (j==0) mu_s(i) += varmat(i, j) * beta_s(ind_loc(i));
      mu_depth(i) += varmat(i, j) * beta_depth(j) * log(depth(i));
    }
    mu(i) = mu_fixed(i) + mu_s(i) + mu_depth(i); 
    // jnll -= dnorm(weight(i), exp(mu(i)), exp(log_phi), true);
    if (!isNA(weight(i))){
      jnll -= dnorm(log(weight(i)), mu(i)-(exp(log_phi)*exp(log_phi))/2, exp(log_phi),true);
    }
    // jnll -= dgamma(weight(i), exp(mu(i))*exp(mu(i))/exp(log_phi), exp(log_phi)/exp(mu(i)), true);
  }
  
  SIMULATE{
    for (int i = 0; i < nobs; i++){
      weight(i) = exp(rnorm(mu(i)-(exp(log_phi)*exp(log_phi))/2,exp(log_phi)));
    }
    REPORT(weight);
  }
  
  
  // Report
  REPORT(beta);
  REPORT(beta_s);
  REPORT(beta_b_s);
  REPORT(beta_depth);
  REPORT(mu);
  REPORT(mu_fixed);
  REPORT(mu_s);
  REPORT(mu_depth);
  ADREPORT(mu);
  ADREPORT(beta);
  ADREPORT(beta_s);
  ADREPORT(beta_b_s);
  ADREPORT(beta_depth);
  ADREPORT(log_phi);
  ADREPORT(log_rho);
  ADREPORT(log_sig);
  ADREPORT(log_rho_b);
  ADREPORT(log_sig_b);
  
  
  return(jnll);
  
}
