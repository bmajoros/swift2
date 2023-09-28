/****************************************************************
 Swift.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Swift.H"
#include "BOOM/GSL/BetaDistribution.H"
using namespace std;
using namespace BOOM;


Swift::Swift(float conc)
  : conc(conc)
{
  // ctor
}



void Swift::run(const Replicates &DNA,const Replicates &RNA,
		const int numSamples,Array1D<SwiftSample> &samples)
{
  if(samples.size()!=numSamples) samples.resize(numSamples);
  const int numDnaReps=DNA.size(), numRnaReps=RNA.size();
  int nextDnaRep=0, nextRnaRep=0;
  for(int i=0 ; i<numSamples ; ++i) {
    // Get counts from current reps
    const int a=DNA[nextDnaRep].getAlt(), b=DNA[nextDnaRep].getRef();
    const int k=RNA[nextRnaRep].getAlt(), m=RNA[nextRnaRep].getRef();
    
    // Sample p from P(p|a,b)
    GSL::BetaDistribution beta1(a+1,b+1); // posterior with uniform prior
    float p;
    do { p=beta1.random(); } while(p==0.0 || p==1.0);

    // Sample q from P(q|p,k,m)
    float alpha, beta;
    chooseAlphaBeta(conc,p,alpha,beta);
    GSL::BetaDistribution beta2(k+alpha,m+beta); // posterior
    float q;
    do { q=beta2.random(); } while(q==0.0 || q==1.0);

    // Create a new sample via (q/(1-q))/(p/(1-p))
    samples[i]=SwiftSample(p,q);

    // Advance to next pairing of DNA & RNA replicates
    advanceReps(nextDnaRep,nextRnaRep,numDnaReps,numRnaReps);
  }
}



void Swift::chooseAlphaBeta(float c,const float p,float &alpha,float &beta)
{
  // If c>2, this will use a beta prior with p as the mode.
  // For c=2, this will result in a uniform prior.
  if(c<2) c=2;

  // ### ADAPTIVE CONCENTRATION METHOD:
  //const float c=5.5/p;
  // ###

  // Compute parameters for a beta prior with mode and concentration
  alpha=p*(c-2)+1;
  beta=(1-p)*(c-2)+1;
}



void Swift::advanceReps(int &dnaIndex,int &rnaIndex,const int numDnaReps,
			const int numRnaReps)
{
  ++dnaIndex;
  if(dnaIndex>=numDnaReps) {
    dnaIndex=0;
    rnaIndex=(rnaIndex+1)%numRnaReps;
  }
}




