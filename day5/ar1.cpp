#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER(logSigma);
  PARAMETER(phiTrans);
  PARAMETER_VECTOR(gamma);
  
  Type phi =  Type(2)/(1 + exp(-2*phiTrans))-Type(1);
  Type sd = exp(logSigma); 
  Type nll=0;  
  nll+= -dnorm(gamma(0),Type(0),sqrt(sd*sd/(Type(1)-phi*phi)),true);
  //random effect part
  for(int i =1; i<gamma.size(); ++i){
    nll += -dnorm(gamma(i),phi*gamma(i-1),sd,true);
  }
  // observation part
  for(int i=1;i<y.size();i++){    
    nll += -keep(i)*dpois(y(i),exp(gamma(i)),true);
  }
  REPORT(phi);
  REPORT(sd);
  return nll;
}
