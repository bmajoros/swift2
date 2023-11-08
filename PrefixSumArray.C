/****************************************************************
 PrefixSumArray.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "PrefixSumArray.H"
using namespace std;
using namespace BOOM;


PrefixSumArray::PrefixSumArray(const Trapezoids &trapezoids)
  : array(trapezoids.size()+1)
{
  // ctor

  compute(trapezoids);
}



int PrefixSumArray::size() const
{
  return array.size();
}



const float &PrefixSumArray::operator[](int i) const
{
  return array[i];
}



void PrefixSumArray::compute(const Trapezoids &trapezoids)
{
  const int numTrap=trapezoids.size();
  const int psaSize=array.size();
  array[0]=0.0;
  for(int i=1 ; i<psaSize ; ++i)
    array[i]=array[i-1]+trapezoids[i-1];
}



