/****************************************************************
 InverseCDF.H
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_InverseCDF_H
#define INCL_InverseCDF_H
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Trapezoids.H"
#include "GridMap.H"
using namespace std;
using namespace BOOM;


class InverseCDF {
public:
  InverseCDF(const Trapezoids &,const GridMap &);
  virtual ~InverseCDF();
  float sample();
private:
  const GridMap &map;
  gsl_ran_discrete_t *gsl_table;
  void computeTable(const Trapezoids &);
};

#endif

