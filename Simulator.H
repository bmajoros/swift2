/****************************************************************
 Simulator.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Simulator_H
#define INCL_Simulator_H
#include <iostream>
#include "Experiment.H"
#include "BOOM/Array1D.H"
#include "BOOM/GSL/BetaDistribution.H"
using namespace std;
using namespace BOOM;


class Simulator {
public:
  void sim(const Experiment &realData,const int numNulls,
	   Array1D<Experiment> &nulls);
private:
  Array1D<int> dnaSums, rnaSums;
  Array1D<GSL::BetaDistribution> betas; // for DNA replicates
  void computeSums(const Experiment &realData);
  void computeSums(const Replicates &data,Array1D<int> &into);
  void prepareBetas(const Replicates &DNA);
};

#endif

