/****************************************************************
 DensityGrid.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "DensityGrid.H"
#include "GridMap.H"
using namespace std;
using namespace BOOM;

DensityGrid::DensityGrid(int numPoints)
  : densities(numPoints)
{
  // ctor
}



int DensityGrid::size() const
{
  return densities.size();
}



float &DensityGrid::operator[](int i)
{
  return densities[i];
}



const float &DensityGrid::operator[](int i) const
{
  return densities[i];
}



void DensityGrid::fill(const DensityFunction &f)
{
  const int n=size();
  GridMap gridMap(n);
  for(int i=0 ; i<n ; ++i) {
    const float p=gridMap.indexToP(i);
    if(p==0 || p==1) densities[i]=0;
    else {
      densities[i]=f(p);
      //cout<<"p="<<p<<", f(p)="<<f(p)<<endl;
    }
  }
}



