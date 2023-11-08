/****************************************************************
 DensityGrid.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_DensityGrid_H
#define INCL_DensityGrid_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "DensityFunction.H"
using namespace std;
using namespace BOOM;

class DensityGrid {
public:
  DensityGrid(int numPoints);
  int size() const;
  float &operator[](int);
  const float &operator[](int) const;
  void fill(const DensityFunction &);
private:
  Array1D<float> densities;
};

#endif

