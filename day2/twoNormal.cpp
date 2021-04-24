#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_ARRAY(red);
  DATA_ARRAY(blue);
  DATA_ARRAY(black);

  PARAMETER_VECTOR(muRed);
  PARAMETER_VECTOR(logSigmaRed);
  PARAMETER(logitRhoRed);
  vector<Type> sigmaRed = exp(logSigmaRed);
  Type rhoRed=2*exp(logitRhoRed)/(exp(logitRhoRed)+1)-1;
  matrix<Type> SigmaRed(2,2);
  SigmaRed(0,0)=pow(sigmaRed(0),2); 
  SigmaRed(0,1)=sigmaRed(0)*sigmaRed(1)*rhoRed; 
  SigmaRed(1,0)=SigmaRed(0,1); 
  SigmaRed(1,1)=pow(sigmaRed(1),2); 
  MVNORM_t<Type> nldensRed(SigmaRed);

  PARAMETER_VECTOR(muBlue);
  PARAMETER_VECTOR(logSigmaBlue);
  PARAMETER(logitRhoBlue);
  vector<Type> sigmaBlue = exp(logSigmaBlue);
  Type rhoBlue=2*exp(logitRhoBlue)/(exp(logitRhoBlue)+1)-1;
  matrix<Type> SigmaBlue(2,2);
  SigmaBlue(0,0)=pow(sigmaBlue(0),2); 
  SigmaBlue(0,1)=sigmaBlue(0)*sigmaBlue(1)*rhoBlue; 
  SigmaBlue(1,0)=SigmaBlue(0,1); 
  SigmaBlue(1,1)=pow(sigmaBlue(1),2); 
  MVNORM_t<Type> nldensBlue(SigmaBlue);

  Type nll = 0; 
  for(int i=0; i<red.dim(1); ++i){
    nll+=nldensRed(red.col(i)-muRed);
  }
  for(int i=0; i<blue.dim(1); ++i){
    nll+=nldensBlue(blue.col(i)-muBlue);
  }

  vector<int> res(black.dim(1));  
  for(int i=0; i<black.dim(1); ++i){
    if(nldensRed(black.col(i)-muRed) < nldensBlue(black.col(i)-muBlue) ){
      res(i)=0;
    }else{
      res(i)=1;
    }
  }

  REPORT(res)
  return nll;
}
