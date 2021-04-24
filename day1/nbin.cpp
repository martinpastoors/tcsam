#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  PARAMETER(logsize);
  PARAMETER(p);
  Type size = exp(logsize);
  Type nll = -sum(dnbinom(Y, size, p, true)); // test
  ADREPORT(size);
  return nll;
}
