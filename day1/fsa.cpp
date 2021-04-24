#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_IVECTOR(year)
  DATA_IVECTOR(fleet)
  DATA_IVECTOR(age)
  DATA_VECTOR(obs)
  DATA_ARRAY(stockMeanWeight)
  DATA_ARRAY(M)              
  DATA_ARRAY(propMature)     
  DATA_SCALAR(surveyTime)     

  PARAMETER_VECTOR(logN1Y);   // read log N in the first year
  PARAMETER_VECTOR(logN1A);   // read log N in the first age
  PARAMETER_VECTOR(logFY);    // read log F by year
  PARAMETER_VECTOR(logFA);    // read log F by age
  PARAMETER(logVarLogCatch);  // Log variance of log catch
  PARAMETER_VECTOR(logQ);     // Log of catchability parameters
  PARAMETER(logVarLogSurvey); // Lo variance of the log survey 

  int minAge=age.minCoeff();  // minimum age
  int maxAge=age.maxCoeff();  // maximum age
  int na=maxAge-minAge+1;     // max - min +1
  int minYear=year.minCoeff();// min year etc
  int maxYear=year.maxCoeff();// max year
  int ny=maxYear-minYear+1;   // number of years

  Type ans=0;   // answer

  // setup F matrix
  matrix<Type> F(na,ny);        // set up F matrix
  for(int a=0; a<na; ++a){      // increment age groups ++a = increment by 1
    for(int y=0; y<ny; ++y){
      F(a,y)=exp(logFY(y))*exp(logFA(a));
    }
  }
  
  // setup logN
  matrix<Type> logN(na,ny);
  for(int a=0; a<na; ++a){
    logN(a,0)=logN1Y(a);        // fill the first year; all ages
  } 
  for(int y=1; y<ny; ++y){
    logN(0,y)=logN1A(y-1);
    for(int a=1; a<na; ++a){
      logN(a,y)=logN(a-1,y-1)-F(a-1,y-1)-M(a-1,y-1);
    }
  } 

  // Match to observations
  vector<Type> logObs=log(obs);  // log of all observations
  Type pred, sd;
  int a, y;
  for(int i=0; i<logObs.size(); ++i){
    a = age(i)-minAge;           // calculate age
    y = year(i)-minYear;         // calculate year

    if(fleet(i)==1){
      pred = log(F(a,y))-log(F(a,y)+M(a,y))+log(Type(1.0)-exp(-F(a,y)-M(a,y)))+logN(a,y);
      sd = exp(Type(0.5)*logVarLogCatch);
    }else{
      pred = logQ(a)-(F(a,y)+M(a,y))*surveyTime+logN(a,y);
      sd = exp(Type(0.5)*logVarLogSurvey);
    }    
    ans += -dnorm(logObs(i),pred,sd,true);  // add one at a time
  }

  vector<Type> ssb(ny);
  ssb.setZero();
  for(int y=0; y<ny; ++y){
    for(int a=0; a<na; ++a){
      ssb(y)+=exp(logN(a,y))*stockMeanWeight(a,y)*propMature(a,y);
    }
  }

  ADREPORT(ssb);
  ADREPORT(F);
  return ans;
}
