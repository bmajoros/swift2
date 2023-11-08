/****************************************************************
 PosteriorEstimator.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "PosteriorEstimator.H"
using namespace std;
using namespace BOOM;


PosteriorEstimator::PosteriorEstimator(const PrefixSumArray &psa)
  : psa(psa)
{
  //ctor
}



float PosteriorEstimator::operator()(int i,int j) const
{
  const int n=psa.size()-1;
  const float total=psa[n];
  const float numer=psa[i]+total-psa[j];
  const float posterior=numer/total;
  return posterior;
}


float PosteriorEstimator::totalArea() const
{
  return psa[psa.size()-1];
}




