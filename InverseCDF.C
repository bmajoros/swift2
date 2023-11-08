/****************************************************************
 InverseCDF.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "InverseCDF.H"
#include "BOOM/GSL/Random.H"
using namespace std;
using namespace BOOM;
using namespace GSL;

InverseCDF::InverseCDF(const Trapezoids &trapezoids,const GridMap &map)
  : map(map)
{
  computeTable(trapezoids);
}



InverseCDF::~InverseCDF()
{
  gsl_ran_discrete_free(gsl_table);
}



void InverseCDF::computeTable(const Trapezoids &trapezoids)
{
  // First, compute total area
  const int N=trapezoids.size();
  float total=0;
  for(int i=0 ; i<N ; ++i) total+=trapezoids[i];
  
  // Now normalize trapezoids by total area
  Array1D<double> normalized(N);
  for(int i=0 ; i<N ; ++i) normalized[i]=trapezoids[i]/total;

  // Give normalized probs to GSL to make its table
  gsl_table=gsl_ran_discrete_preproc(N,normalized.getRawArray());
}



float InverseCDF::sample()
{
  const gsl_rng *rng=Random::getGenerator();
  size_t index=gsl_ran_discrete(rng,gsl_table);
  const float p=map.indexToP(index);
  return p;
}





