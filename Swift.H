/****************************************************************
 Swift.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Swift_H
#define INCL_Swift_H
#include <iostream>
using namespace std;
using namespace BOOM;


class Swift {
public:
  Swift(const float conc);
  void run(const Replicates &DNA,const Replicates &RNA,
	   Array1D<float> &samples);
};

#endif

