/****************************************************************
 PrefixSumArray.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_PrefixSumArray_H
#define INCL_PrefixSumArray_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "Trapezoids.H"
using namespace std;
using namespace BOOM;


class PrefixSumArray {
public:
  PrefixSumArray(const Trapezoids &);
  int size() const;
  const float &operator[](int) const;
private:
  Array1D<float> array;
  void compute(const Trapezoids &);
};

#endif

