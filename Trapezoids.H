/****************************************************************
 Trapezoids.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Trapezoids_H
#define INCL_Trapezoids_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "DensityGrid.H"
using namespace std;
using namespace BOOM;

class Trapezoids {
public:
  Trapezoids(const DensityGrid &);
  int size() const;
  const float &operator[](int) const;
private:
  Array1D<float> areas;
  void computeAreas(const DensityGrid &);
};

#endif

