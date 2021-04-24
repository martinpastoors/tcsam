#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  DATA_VECTOR(y);
  DATA_INTEGER(code);
  PARAMETER(logSigma);
  PARAMETER(phiTrans);
  PARAMETER_VECTOR(gamma);
  
  Type phi =  Type(2)/(1 + exp(-2*phiTrans))-Type(1);
  Type sd = exp(logSigma); 
    
  Type nll=0;
  
  
  if(code==0){
    nll+= -dnorm(gamma(0),Type(0),sqrt(sd*sd/(Type(1)-phi*phi)),true);
    for(int i =1; i<gamma.size(); ++i){
      nll += -dnorm(gamma(i),phi*gamma(i-1),sd,true);
    }
  }
  if(code==1){
    nll+=SCALE(AR1(phi),sqrt(sd*sd/(1-phi*phi)))(gamma);
  }
  
  for(int i=1;i<y.size();i++){    
    nll += -dpois(y(i),exp(gamma(i)),true);
  }
  return nll;
}




