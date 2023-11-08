/****************************************************************
 GridMap.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_GridMap_H
#define INCL_GridMap_H
#include <iostream>
using namespace std;

class GridMap {
public:
  GridMap(int gridSize);
  int pToIndex(float p) const;
  float indexToP(int index) const;
  int getGridSize() const { return gridSize; }
private:
  const int gridSize;
};

#endif

