/****************************************************************
 DensityWithPrior.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityWithPrior_H
#define INCL_DensityWithPrior_H
#include <iostream>
#include "DensityFunction.H"
using namespace std;
using namespace BOOM;

// This is a version of the DensityFunction that includes a shrinage
// prior.  Given sample allele frequency in DNA (p), this shrinks the
// allele frequency in RNA (q) toward p.
class DensityWithPrior : public DensityFunction {
public:
  DensityWithPrior(const Replicates &,const float repConc,const float p);
  virtual double operator()(float q) const; // Not in log space!
  void changeP(const float p) { this->p=p; }
protected:
  float p;
};

#endif

