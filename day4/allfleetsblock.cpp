#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

bool isNAINT(int x){
  return NA_INTEGER==x;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(N)
  DATA_MATRIX(F)
  DATA_MATRIX(M)
  DATA_IMATRIX(aux)
  DATA_IMATRIX(idx1)
  DATA_IMATRIX(idx2)   
  DATA_IVECTOR(fleetTypes)
  DATA_VECTOR(sampleTimes)    
  DATA_INTEGER(minYear)
  DATA_INTEGER(minAge)
  DATA_VECTOR(obs)
  DATA_IMATRIX(keyQ)   
  DATA_IMATRIX(keySd)
  int nobs=obs.size();
  
  PARAMETER_VECTOR(logQ)
  PARAMETER_VECTOR(logsd)
  PARAMETER_VECTOR(missing)
  vector<Type> sd=exp(logsd);

  //patch missing 
  int idxmis=0; 
  for(int i=0; i<nobs; i++){
    if(isNA(obs(i))){
      obs(i)=exp(missing(idxmis++));
    }    
  }

  Type nll=0;
  
  vector<Type> logPred(nobs);
  Type Z;
  int y, a, f;
  for(int i=0; i<nobs; ++i){
    y=aux(i,0)-minYear;
    f=aux(i,1)-1;
    a=aux(i,2)-minAge;
    Z=F(y,a)+M(y,a);
    switch(fleetTypes(f)){
      case 0:
        logPred(i) = log(N(y,a))-log(Z)+log(1-exp(-Z))+log(F(y,a));
      break;
      
      case 2:
        logPred(i) = logQ(keyQ(f,a))+log(N(y,a))-Z*sampleTimes(f);
      break;

      default:
        std::cout<<"Error: This fleet type not implemented yet."<<std::endl;
        exit(EXIT_FAILURE);
      break;
    }
  }

  vector< density::MVNORM_t<Type> >  nllVec(fleetTypes.size());

  for(int f=0; f<idx1.rows(); ++f){
    int d=0;
    for(int a=0; a<keySd.cols(); ++a){
      if(!isNAINT(keySd(f,a)))d++;
    };
    matrix<Type> S(d,d);
    S.setZero();
    d=0;
    for(int a=0; a<keySd.cols(); ++a){
      if(!isNAINT(keySd(f,a))){
	S(d,d)=sd(keySd(f,a))*sd(keySd(f,a));
	d++;
      }
    };
    nllVec(f).setSigma(S);
  }
  
  for(int f=0; f<idx1.rows(); ++f){
    for(int y=0; y<idx1.cols(); ++y){
      if(!isNAINT(idx1(f,y))){
        vector<Type> o=obs.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
	vector<Type> p=logPred.segment(idx1(f,y),idx2(f,y)-idx1(f,y)+1);
	nll += nllVec(f)(log(o)-p);
      }
    }
  }

  REPORT(logPred)
  return nll;
}
