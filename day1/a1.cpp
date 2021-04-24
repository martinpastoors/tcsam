#include<TMB.hpp>

template <class Type>
Type objective_function <Type> :: operator() ()
{
  PARAMETER(x);
  return (x-Type(42.0))*(x-Type(42.0));
}

