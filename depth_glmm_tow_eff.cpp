// ------------------------------------------------------------------------
// modeling of weight~height+age+... : spatial slope
// ------------------------------------------------------------------------

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() (){
  
  using namespace density;
  
  
  DATA_VECTOR(weight);
  DATA_MATRIX(varmat);
  DATA_VECTOR(depth);
  DATA_IVECTOR(tow_id);
  
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(beta_depth);
  PARAMETER(log_phi);
  PARAMETER(log_epsilon);
  
  PARAMETER_VECTOR(tow_eff);
  
  // ---------------------------
  // Joint negative log-likelihood
  Type jnll = Type(0.0);
  
  int n_tow = tow_eff.size();
  for (int tow = 0; tow < n_tow; tow++){
    jnll -= dnorm(tow_eff(tow), Type(0.0) ,exp(log_epsilon),true);
  }
  
  // Spatial effect: Random intercept and Random slope
  
  // Observation likelihood
  int nobs = weight.size();
  vector<Type> mu(nobs); mu.setZero(); 
  vector<Type> mu_fixed(nobs); mu_fixed.setZero(); 
  vector<Type> mu_depth(nobs); mu_depth.setZero();  
  vector<Type> mu_tow(nobs); mu_tow.setZero();
  for (int i = 0; i < nobs; i++){
    for(int j = 0; j < varmat.cols(); j++){
      mu_fixed(i) += varmat(i, j) * beta(j);
      mu_depth(i) += varmat(i, j) * beta_depth(j) * log(depth(i));
      if (j==0) mu_tow(i) += tow_eff(tow_id(i));
    }
    mu(i) = mu_fixed(i) + mu_depth(i) + mu_tow(i); 
    // jnll -= dnorm(weight(i), exp(mu(i)), exp(log_phi), true);
    jnll -= dnorm(log(weight(i)), mu(i)-(exp(log_phi)*exp(log_phi))/2, exp(log_phi),true);
    // jnll -= dgamma(weight(i), exp(mu(i))*exp(mu(i))/exp(log_phi), exp(log_phi)/exp(mu(i)), true);
  }
  
  
  
  // Report
  REPORT(beta);
  REPORT(beta_depth);
  REPORT(mu);
  REPORT(mu_fixed);
  REPORT(mu_depth);
  REPORT(mu_tow);
  // ADREPORT(mu);
  ADREPORT(beta);
  ADREPORT(beta_depth);
  ADREPORT(log_phi);
  
  
  return(jnll);
  
}
