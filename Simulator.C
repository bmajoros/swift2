/****************************************************************
 Simulator.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Simulator.H"
#include "BOOM/Vector.H"
#include "BOOM/GSL/GslBinomial.H"
using namespace std;
using namespace BOOM;


void Simulator::sim(const Experiment &realData,const int numNulls,
		    Array1D<Experiment> &nulls)
{
  // Prepare output and precompute some things
  const int numDnaReps=realData.DNA.size(), numRnaReps=realData.RNA.size();
  computeSums(realData);
  prepareBetas(realData.DNA);
  nulls.resize(numNulls);
  for(int i=0 ; i<numNulls ; ++i) {
    Experiment &null=nulls[i];
    null.DNA.resize(numDnaReps);
    null.RNA.resize(numRnaReps);
  }

  // Perform simulation
  GSL::GslBinomial binom;
  GSL::BetaDistribution beta;
  for(int i=0 ; i<numNulls ; ++i) {
    // First, sample a value for p from the real DNA (using a beta distr)
    // and use it to simulate new DNA allele counts
    Experiment &null=nulls[i];
    Vector<float> Ps;
    for(int j=0 ; j<numDnaReps ; ++j) {
      const float p=/*0.5;*/betas[j].random(); // ### DEBUGGING: p=0.5
      Ps.push_back(p);
      binom.change(p);
      const int n=dnaSums[j];
      const int a=binom.random(n);
      const int b=n-a;
      Replicate &dnaRep=null.DNA[j];
      dnaRep.setRef(a); dnaRep.setAlt(b);
    }
    
    // Assuming q=p, simulate RNA counts from a binomial
    int nextP=0;
    for(int j=0 ; j<numRnaReps ; ++j) {
      //const float q=Ps[nextP];

      // ### EXPERIMENTAL:
      const float mode=Ps[nextP];
      const float conc=139;
      const float ALPHA=mode*(conc-2)+1;
      const float BETA=(1-mode)*(conc-2)+1;
      beta.change(ALPHA,BETA);
      const float q=beta.random();
      // ###
      
      nextP=(nextP+1)%Ps.size();
      binom.change(q);
      const int n=rnaSums[j];
      const int k=binom.random(n);
      const int m=n-k;
      Replicate &rnaRep=null.RNA[j];
      rnaRep.setRef(k); rnaRep.setAlt(m);
    }
  }
}



void Simulator::prepareBetas(const Replicates &DNA)
{
  const int numReps=DNA.size();
  betas.resize(numReps);
  for(int i=0 ; i<numReps ; ++i)
    betas[i].change(DNA[i].getRef()+1,DNA[i].getAlt()+1);
}



// Array1D<int> dnaSums, rnaSums;
void Simulator::computeSums(const Experiment &realData)
{
  computeSums(realData.DNA,dnaSums);
  computeSums(realData.RNA,rnaSums);
}



void Simulator::computeSums(const Replicates &data,Array1D<int> &into)
{
  const int numReps=data.size();
  into.resize(numReps);
  for(int i=0 ; i<numReps ; ++i) {
    into[i]=data[i].getSum();
  }
}





  
