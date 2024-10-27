// ------------------------------------------------------------------------
// modeling of weight~height+age+... : spatial slope
// ------------------------------------------------------------------------

#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() (){
  
  using namespace density;
  
  
  DATA_VECTOR(weight);
  DATA_VECTOR(heights);
  DATA_FACTOR(tow_id);
  
  PARAMETER(b);
  PARAMETER(beta);
  PARAMETER(log_phi);
  PARAMETER(log_epsilon);
  
  PARAMETER_VECTOR(tow_eff);
  
  // ---------------------------
  // Joint negative log-likelihood
  Type jnll = Type(0.0);
  
  // Tow effect
  int n_tow = tow_eff.size();
  for (int tow = 0; tow < n_tow; tow++){
    jnll -= dnorm(tow_eff(tow), Type(0.0) ,exp(log_epsilon),true);
  }
  
  // Observation likelihood
  int nobs = weight.size();
  vector<Type> mu(nobs); mu.setZero(); 
  for (int i = 0; i < nobs; i++){
    mu(i) =  (beta + tow_eff(tow_id(i))) * exp(b*log(heights(i)));
    jnll -= dnorm(weight(i), mu(i), exp(log_phi), true);
    // jnll -= dnorm(log(weight(i)), mu(i)-(exp(log_phi)*exp(log_phi))/2, exp(log_phi),true);
    // jnll -= dgamma(weight(i), exp(mu(i))*exp(mu(i))/exp(log_phi), exp(log_phi)/exp(mu(i)), true);
  }
  
  SIMULATE{
    for (int i = 0; i < nobs; i++){
      weight(i) = rnorm(mu(i),exp(log_phi));
    }
  }
  
  
  // Report
  REPORT(beta);
  REPORT(mu);
  REPORT(tow_eff);
  // ADREPORT(mu);
  ADREPORT(b);
  ADREPORT(beta);
  ADREPORT(log_phi);
  ADREPORT(log_epsilon);
  ADREPORT(tow_eff);
  
  
  return(jnll);
  
}
