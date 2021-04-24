#include <TMB.hpp>
 
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(SSB)
  DATA_VECTOR(logR);

  PARAMETER(logA); 
  PARAMETER(logB);
  PARAMETER(logSigma);
  vector<Type> pred = logA + log(SSB) - log(Type(1)+exp(logB) * SSB);
  Type ans=-sum(dnorm(logR, pred, exp(logSigma), true));
  ADREPORT(pred)
  return ans;
}
