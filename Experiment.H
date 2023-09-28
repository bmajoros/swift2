/****************************************************************
 Experiment.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Experiment_H
#define INCL_Experiment_H
#include <iostream>
#include "Replicates.H"
using namespace std;
using namespace BOOM;

struct Experiment {
public:
  Replicates DNA, RNA;
};

#endif

