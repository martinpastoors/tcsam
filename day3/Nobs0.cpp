#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(year)
  DATA_VECTOR(age)
  DATA_MATRIX(M)
  DATA_MATRIX(F)
  DATA_MATRIX(Nobs)
  int mode=0;  
  int nrow=Nobs.rows();
  int ncol=Nobs.cols();
  
  PARAMETER(logsdR)
  PARAMETER(logsdS)
  PARAMETER(logsd)
  PARAMETER_MATRIX(logN)
  Type sdR = exp(logsdR);
  Type sdS = exp(logsdS);
  Type sd = exp(logsd);  

  Type jnll=0;
  
  Type pred;
  for(int y=1; y<nrow; ++y){
    switch(mode){
      case 0: // RW
        pred = logN(y-1,0);
      break;

      case 1:
	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
      break;

      case 2:
	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
      break;

      default:
	std::cout<<"Stock recruitment code not implemented yet."<<std::endl;
      break;
    }
    jnll += -dnorm(logN(y,0),pred,sdR,true);
  }
  
  // TO DO: MAKE SURE THAT THE Y,A SEQUENCE IS USED; or invert the data before it enters into C
  
  for(int y=1; y<nrow; ++y){
    for(int a=1; a<ncol; ++a){
      pred=logN(y-1,a-1)-F(y-1,a-1)-M(y-1,a-1);
      if(a==(ncol-1)){
        pred=log(exp(pred) +  exp(logN(y-1,a)-F(y-1,a)-M(y-1,a)));
      }
      jnll += -dnorm(logN(y,a),pred,sdS,true);
    }
  } 

  for(int y=0; y<nrow; ++y){
    for(int a=0; a<ncol; ++a){
      jnll += -dnorm(Nobs(y,a), logN(y,a),sd,true);
    }
  } 
  
  return jnll;
}


