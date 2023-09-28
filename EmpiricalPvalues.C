/****************************************************************
 EmpiricalPvalues.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "EmpiricalPvalues.H"
using namespace std;
using namespace BOOM;


float EmpiricalPvalues::median_p(const Array1D<SwiftSample> &realSamples,
			 const Array1D< Array1D<SwiftSample> > &nullSamples)
{
}



float EmpiricalPvalues::area_p(const Array1D<SwiftSample> &realSamples,
	     const float lambda, // e.g., 1.2
	     const Array1D< Array1D<SwiftSample> > &nullSamples)
{
}



float EmpiricalPvalues::getMedian(const Array1D<SwiftSample> &)
{
}



float EmpiricalPvalues::getArea(const Array1D<SwiftSample> &,
				const float lambda)
{
}









