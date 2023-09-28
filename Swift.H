/****************************************************************
 Swift.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Swift_H
#define INCL_Swift_H
#include <iostream>
#include "Replicates.H"
#include "BOOM/Array1D.H"
#include "SwiftSample.H"
using namespace std;
using namespace BOOM;


class Swift {
public:
  Swift(const float conc);
  void run(const Replicates &DNA,const Replicates &RNA,const int numSamples,
	   Array1D<SwiftSample> &samples);
private:
  float conc;
  void chooseAlphaBeta(float c,const float p,float &alpha,float &beta);
  void advanceReps(int &dnaIndex,int &rnaIndex,const int numDnaReps,
		   const int numRnaReps);
};

#endif

