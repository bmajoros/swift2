/****************************************************************
 PosteriorEstimator.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_PosteriorEstimator_H
#define INCL_PosteriorEstimator_H
#include <iostream>
#include "PrefixSumArray.H"
using namespace std;
using namespace BOOM;

class PosteriorEstimator {
public:
  PosteriorEstimator(const PrefixSumArray &);
  float operator()(int i,int j) const; // indices into the GridArray
  float totalArea() const;
private:
  const PrefixSumArray &psa;
};

#endif

