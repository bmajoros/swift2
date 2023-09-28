/****************************************************************
 EmpiricalPvalues.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_EmpiricalPvalues_H
#define INCL_EmpiricalPvalues_H
#include <iostream>
#include "SwiftSample.H"
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;

class EmpiricalPvalues {
public:
  float median_p(Array1D<SwiftSample> &realSamples,
		 Array1D< Array1D<SwiftSample> > &nullSamples);
  float area_p(Array1D<SwiftSample> &realSamples,
	       const float lambda, // e.g., 1.2
	       Array1D< Array1D<SwiftSample> > &nullSamples);
private:
  bool isSorted(const Array1D<SwiftSample> &);
  void sort(Array1D<SwiftSample> &);
  float getMedian(Array1D<SwiftSample> &);
  float getArea(Array1D<SwiftSample> &,const float lambda);
};

#endif

