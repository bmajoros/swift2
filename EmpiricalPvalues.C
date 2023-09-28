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


bool EmpircalPvalues::isSorted(const Array1D<SwiftSample> &samples)
{
  const int n=samples.size();
  for(int i=0 ; i<n-1 ; ++i)
    if(samples[i].getTheta()>samples[i+1].getTheta()) return false;
  return true;
}



void EmpircalPvalues::sort(Array1D<SwiftSample> &samples)
{
  SwiftSampleComparator cmp;
  Array1DSorter<SwiftSample> sorter(samples,cmp);
  sorter.sortAscendInPlace();
}



float EmpiricalPvalues::median_p(Array1D<SwiftSample> &realSamples,
			Array1D< Array1D<SwiftSample> > &nullSamples)
{
  const float realMedian=getMedian(realSamples);

  
}



float EmpiricalPvalues::area_p(Array1D<SwiftSample> &realSamples,
	            const float lambda, // e.g., 1.2
	            Array1D< Array1D<SwiftSample> > &nullSamples)
{
  const float realArea=getArea(realSamples,lambda);


  
}



float EmpiricalPvalues::getMedian(Array1D<SwiftSample> &samples)
{
  if(!isSorted(samples)) sort(samples);
  int n=samples.size();
  if(n<2) throw "Too few samples to identify median";
  int mid=n/2;
  float median;
  if(n%2==0) median=(samples[mid-1].getTheta()+samples[mid].getTheta())/2.0;
  else median=samples[mid].getTheta();
  return median;
}



float EmpiricalPvalues::getArea(Array1D<SwiftSample> &samples,
				const float lambda)
{
  if(!isSorted(samples)) sort(samples);
  const float invLambda=1.0/lambda;
  float count=0;
  const int numSamples=samples.size();
  for(int i=0 ; i<numSamples ; ++i) {
    const theta=samples[i].getTheta();
    if(theta<=invLambda || theta>=lambda) ++count;
  }
  const float p=float(count)/float(numSamples);
  return p;
}









