/****************************************************************
 DensityFunction.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityFunction_H
#define INCL_DensityFunction_H
#include <iostream>
#include <gsl/gsl_sf_gamma.h>
#include "Replicates.H"
using namespace std;
using namespace BOOM;


class DensityFunction {
public:
  DensityFunction(const Replicates &,const float concentration);
  virtual double operator()(float q) const; // Not in log space!
protected:
  const float c; // concentration parm for between-replicate variation
  const Replicates &data;
  inline double lg(double x) { return gsl_sf_lngamma(x); }
};

#endif

