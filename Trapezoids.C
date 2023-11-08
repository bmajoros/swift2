/****************************************************************
 Trapezoids.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Trapezoids.H"
using namespace std;
using namespace BOOM;

Trapezoids::Trapezoids(const DensityGrid &grid)
  : areas(grid.size()-1)
{
  // ctor

  computeAreas(grid);
}



int Trapezoids::size() const
{
  return areas.size();
}



const float &Trapezoids::operator[](int i) const
{
  return areas[i];
}



void Trapezoids::computeAreas(const DensityGrid &grid)
{
  const int numTrapezoids=areas.size();
  const float dx=1.0/float(numTrapezoids);
  for(int i=0 ; i<numTrapezoids ; ++i)
    areas[i]=dx*(grid[i]+grid[i+1])/2.0;
}




