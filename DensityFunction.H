/****************************************************************
 DensityFunction.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityFunction_H
#define INCL_DensityFunction_H
#include <iostream>
#include "Replicates.H"
using namespace std;
using namespace BOOM;


class DensityFunction {
public:
  DensityFunction(const Replicates &,const float concentration);
  double operator()(float q) const; // Not in log space!
private:
  const float c; // concentration
  const Replicates &data;
};

#endif

