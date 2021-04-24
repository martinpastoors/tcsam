#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(N)
  DATA_MATRIX(F)
  DATA_MATRIX(M)
  DATA_IMATRIX(aux)
  DATA_INTEGER(minYear)
  DATA_INTEGER(minAge)
  DATA_VECTOR(Cobs)
    
  int nobs=Cobs.size();
  PARAMETER_VECTOR(logAlpha);
  PARAMETER_VECTOR(logBeta);
  
  vector<Type> alpha=exp(logAlpha);
  vector<Type> beta = exp(logBeta);
  
  Type nll=0;
  
  vector<Type> logPred(nobs);
  Type Z;
  int y, a;
  for(int i=0; i<nobs; ++i){
    y=aux(i,0)-minYear;
    a=aux(i,2)-minAge;
    Z=F(y,a)+M(y,a);
    logPred(i) = log(N(y,a))-log(Z)+log(1-exp(-Z))+log(F(y,a));
    Type sdThis = sqrt(log (alpha(a)*pow(exp(logPred(i)), (beta(a)-Type(2))) + Type(1)) );
    nll += -dnorm(log(Cobs(i)),logPred(i),sdThis ,true);
  }
  REPORT(logPred)
  return nll;
}
