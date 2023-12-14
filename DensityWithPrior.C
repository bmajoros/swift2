/****************************************************************
 DensityWithPrior.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "DensityWithPrior.H"
using namespace std;
using namespace BOOM;


DensityWithPrior::DensityWithPrior(const Replicates &reps,
				   const float repConc,const float p)
  : DensityFunction(reps,repConc), p(p)
{
  // ctor
}



double DensityWithPrior::operator()(float q) const
{
  // Returned value is not in log space!

  // These are the parameters of the prior on overall freq: P(q|p)
  const float c=5.5/p; // the "adaptive concentration method" from Swift v1
  const float alphaQ=p*(c-2);
  const float betaQ=(1-p)*(c-2);

  // This is the prior P(q|p):
  double logSum=lg(alphaQ+betaQ) - lg(alphaQ) - lg(betaQ)
    +(alphaQ-1)*log(q)+(betaQ-1)*log(1-q);
  
  // These are the parameters of the prior on rep-specific freq: P(qi|q):
  const double alphaI=q*(c-2)+1, betaI=(1-q)*(c-2)+1;

  const int N=data.size(); // number of replicates
  for(int i=0 ; i<N ; ++i) {
    const Replicate &rep=data[i];
    const double ki=rep.getAlt(), mi=rep.getRef();
    const double ni=ki+mi;
    logSum+=lg(alphaI+betaI) +lg(ni+1) +lg(ki+alphaI) +lg(mi+betaI)
      -lg(alphaI) -lg(betaI) -lg(ki+1) -lg(mi+1) -lg(ni+alphaI+betaI);
  }
  return exp(logSum);
}




