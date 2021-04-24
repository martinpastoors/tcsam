#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(year);
  DATA_VECTOR(ssb);
  DATA_VECTOR(Robs);
  DATA_INTEGER(minAge);    
  DATA_INTEGER(mode);
    
  PARAMETER(logsdo);
  PARAMETER(logsdp);
  PARAMETER_VECTOR(logR);
  PARAMETER_VECTOR(rickerpar);
  PARAMETER_VECTOR(bhpar);
  Type sdo = exp(logsdo);
  Type sdp = exp(logsdp);

  Type jnll=0;  // joint negative likelhood (of random effects and observations)

  Type pred;
  for(int i=1; i<logR.size(); ++i){
    switch(mode){   // configuration step
      case 0:       // RW Random Walk
        pred = logR(i-1);
      break;

      case 1: // Ricker
  	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
        pred = rickerpar(0) + log(ssb(i-minAge)) - exp(rickerpar(1))*ssb(i-minAge);
        break;

      case 2: // B-H
  	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
        pred = bhpar(0) + log(ssb(i-minAge)) -log(Type(1)+exp(bhpar(1))*ssb(i-minAge));
        break;

      default:
	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
      break;
    }
    // your code to calculate process likelihood; sum dnorm
    jnll += -dnorm(logR(i), pred, sdp, true);
  }
  for(int i=0; i<Robs.size(); ++i){
    // your code to calculate observation likelihood; distribution of increments 
    jnll += -dnorm(log(Robs(i)), logR(i), sdo, true);
  }
  return jnll;
}
