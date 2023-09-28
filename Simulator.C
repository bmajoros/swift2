/****************************************************************
 Simulator.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Simulator.H"
using namespace std;
using namespace BOOM;


void Simulator::sim(const Experiment &realData,const int numNulls,
		    Array1D<Experiment> &nulls)
{
  computeSums(realData);
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





  
