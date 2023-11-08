/****************************************************************
 DensityFunction.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "DensityFunction.H"
using namespace std;
using namespace BOOM;


inline double lg(double x)
{
  return gsl_sf_lngamma(x);
}



DensityFunction::DensityFunction(const Replicates &data,
				 const float concentration)
  : c(concentration), data(data)
{
  // ctor
}



double DensityFunction::operator()(float q) const
{
  // Returned value is not in log space!

  const double alpha=q*(c-2)+1, beta=(1-q)*(c-2)+1;
  const int N=data.size(); // number of replicates
  double logSum=0;
  for(int i=0 ; i<N ; ++i) {
    const Replicate &rep=data[i];
    const double ki=rep.getAlt(), mi=rep.getRef();
    const double ni=ki+mi;
    logSum+=lg(alpha+beta) +lg(ni+1) +lg(ki+alpha) +lg(mi+beta)
      -lg(alpha) -lg(beta) -lg(ki+1) -lg(mi+1) -lg(ni+alpha+beta);
  }
  return exp(logSum);
}



