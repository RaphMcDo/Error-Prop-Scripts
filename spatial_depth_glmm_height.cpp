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
  
  DATA_MATRIX(varmat_sh);
  DATA_VECTOR(depth_sh);
  DATA_FACTOR(ind_loc_sh);
  DATA_MATRIX(locations);
  DATA_VECTOR(s_heights);
  DATA_VECTOR(a_b_truncation);
  
  //Indicator only to pick out unique means
  DATA_FACTOR(unique_means);
  DATA_INTEGER(n_mod_loc);
  
  //Spatial params
  PARAMETER(log_rho);
  PARAMETER(log_sig);
  PARAMETER(log_nu);
  
  Type rho;rho = exp(log_rho);
  Type sig; sig = exp(log_sig);
  Type nu; nu = exp(log_nu);
  
  PARAMETER_VECTOR(beta_sh_s);
  PARAMETER_VECTOR(beta_depth);
  PARAMETER_VECTOR(beta_sh);
  PARAMETER(log_upsilon);
  
  
  int n_loc = beta_sh_s.size(); 
  int nobs_sh = s_heights.size();
  
  // ---------------------------
  // Joint negative log-likelihood
  Type jnll = Type(0.0);
  
  // Spatial effect: Random intercept and Random slope
  matrix <Type> D(n_loc,n_loc);
  D = DistMat2d(locations);
  matrix<Type> St_sh(n_loc,n_loc);
  St_sh.setZero();
  
  for(int i=0; i<n_loc; ++i){
    St_sh(i,i) = sig*sig;
  }
  for(int i=0; i<n_loc; ++i){
    for(int j=i+1; j<n_loc; ++j){
      St_sh(i,j) = sig*sig*matern(D(i,j), rho, nu); //Matern Covariance Function
      St_sh(j,i) = St_sh(i,j);
    }
  }
  
  jnll += density::MVNORM(St_sh)(beta_sh_s); 
  
  // Observation likelihood
  vector <Type> mu_sh(nobs_sh); mu_sh.setZero(); 
  vector<Type> mu_fixed_sh(nobs_sh); mu_fixed_sh.setZero(); 
  vector <Type> mu_sh_s(nobs_sh); mu_sh_s.setZero();
  vector<Type> mu_depth_sh(nobs_sh); mu_depth_sh.setZero();
  
  vector <Type> residuals(nobs_sh); residuals.setZero();
  
  //For commercial height, it is actually a truncated normal distribution, with truncation at a=100, b= max observation
  
  vector <Type> cdf_norm_div(nobs_sh);
  
  
  //height
  for (int i = 0; i < nobs_sh; i++){
    for(int j = 0; j < varmat_sh.cols(); j++){
      mu_fixed_sh(i) += varmat_sh(i, j) * beta_sh(j);
      if (j==0) mu_sh_s(i) += varmat_sh(i, j) * beta_sh_s(ind_loc_sh(i));
      mu_depth_sh(i) += varmat_sh(i, j) * beta_depth(j) * depth_sh(i);
    }
    mu_sh(i) = mu_fixed_sh(i) + mu_sh_s(i) + mu_depth_sh(i);
    
    cdf_norm_div(i) = pnorm(a_b_truncation(1), mu_sh(i),exp(log_upsilon)) - pnorm(a_b_truncation(0), mu_sh(i),exp(log_upsilon));
    
    if(!isNA(s_heights(i))){
      jnll -= (dnorm(s_heights(i), mu_sh(i), exp(log_upsilon), true)/cdf_norm_div(i)); 
    }
    // jnll -= dnorm(s_heights(i), mu_sh(i), exp(log_upsilon), true);
    // jnll -= dgamma(s_heights(i), exp(mu_sh(i))*exp(log_upsilon), exp(mu_sh(i))*(exp(log_upsilon)*exp(log_upsilon)), true);
    // jnll -= dnorm(log(s_heights(i)), mu_sh_s(i)-(exp(log_upsilon)*exp(log_upsilon))/2, exp(log_upsilon),true);
    // jnll -= dgamma(s_heights(i), exp(mu.col(0)(i))*exp(mu.col(0)(i))/exp(log_upsilon), exp(log_upsilon)/exp(mu.col(0)(i)), true);
    
    residuals(i) = qnorm((pnorm(s_heights(i),mu_sh(i),exp(log_upsilon))-pnorm(a_b_truncation(0),mu_sh(i),exp(log_upsilon)))/cdf_norm_div(i));
  }
  
  SIMULATE{
    for (int i = 0; i < nobs_sh; i++){
      Type temp = 0;
      while (temp < a_b_truncation(0) || temp > a_b_truncation(1)){
        temp = rnorm(mu_sh(i),exp(log_upsilon));
      }
      s_heights(i) = temp;
    }
    REPORT(s_heights);
  }
  
  REPORT(residuals);
  
  //Picking out the unique means based on the unique tows we give to the model
  vector <Type> unique_mu(n_mod_loc);
  for (int loc = 0; loc < n_mod_loc; loc++){
    unique_mu(loc) = mu_sh(unique_means(loc));
  }
  
  // Report
  REPORT(mu_sh);
  REPORT(beta_sh_s);
  REPORT(beta_depth);
  REPORT(beta_sh);
  REPORT(unique_mu);
  ADREPORT(unique_mu);
  ADREPORT(beta_sh_s);
  ADREPORT(beta_depth);
  ADREPORT(beta_sh);
  ADREPORT(log_upsilon);
  ADREPORT(log_rho);
  ADREPORT(log_sig);
  
  
  return(jnll);
  
}
