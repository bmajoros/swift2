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
  float median_p(const Array1D<SwiftSample> &realSamples,
		 const Array1D< Array1D<SwiftSample> > &nullSamples);
  float area_p(const Array1D<SwiftSample> &realSamples,
	       const float lambda, // e.g., 1.2
	       const Array1D< Array1D<SwiftSample> > &nullSamples);
private:
  float getMedian(const Array1D<SwiftSample> &);
  float getArea(const Array1D<SwiftSample> &,const float lambda);
};

#endif

