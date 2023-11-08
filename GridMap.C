/****************************************************************
 GridMap.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "GridMap.H"
using namespace std;


GridMap::GridMap(int gridSize)
  : gridSize(gridSize)
{
  // ctor
}



int GridMap::pToIndex(float p) const
{
  return int(p*gridSize);
}



float GridMap::indexToP(int index) const
{
  return float(index)/float(gridSize-1);
}





