#include <TMB.hpp>

template <class Type>
Type trans(Type x){
  return Type(2)*invlogit(x)-Type(1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  DATA_VECTOR(x);
  DATA_INTEGER(code);
  PARAMETER(logSigma);
  PARAMETER(itPhi);
  Type phi=trans(itPhi);
  int timeSteps=x.size();
  Type sd=exp(logSigma);

  Type ans=0;

  if(code==0){
    ans += -dnorm(x(0),Type(0),sqrt(sd*sd/(1-phi*phi)),true);
    for(int i=1;i<timeSteps;i++){    
      ans += -dnorm(x(i),phi*x(i-1),sd,true);
    }
  }

  if(code==1){
    ans+=SCALE(AR1(phi),sqrt(sd*sd/(1-phi*phi)))(x);
  }

  if(code==2){
    matrix<Type> cov(timeSteps,timeSteps);
    for(int i=0; i<timeSteps; ++i){
      for(int j=0; j<timeSteps; ++j){
        cov(i,j)=sd*sd/(1-phi*phi)*pow(phi,abs(i-j));
      }
    }
    MVNORM_t<Type> nldens(cov);
    ans+=nldens(x);
  }


  if(code==3){
    Eigen::SparseMatrix<Type> Q(timeSteps,timeSteps);
    Q.coeffRef(0,0)=1;
    Q.coeffRef(timeSteps-1,timeSteps-1)=1;
    for(int i=1; i<(timeSteps-1); ++i)Q.coeffRef(i,i)=1+phi*phi;
    for(int i=1; i<timeSteps; ++i){
      Q.coeffRef(i,i-1)=-phi;
      Q.coeffRef(i-1,i)=-phi;
    }
    Q*=Type(1.0)/sd/sd;
    density::GMRF_t<Type> nldens(Q);
    ans+=nldens(x);
  }

  
  ADREPORT(phi); 
  return ans;
}
