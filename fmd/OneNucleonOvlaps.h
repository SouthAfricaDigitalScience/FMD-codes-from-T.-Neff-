/**

   \file OneNucleonOvlaps.h

   calculate spectroscopic amplitudes in momentum space


   (c) 2005-2007 Thomas Neff, Benjamin Hellwig

*/


#ifndef _ONENUCLEONOVLAPS_H
#define _ONENUCLEONOVLAPS_H


#include <stdio.h>

#include "SlaterDet.h"
#include "Projection.h"


typedef struct {
  double qmax;
  int npoints;
  int nalpha;
  int ncosb;
  int recoil;
} OneNucleonOvlapsPara;


#define NSPEC 5
#define DIMSPEC 18

extern ManyBodyOperator OpOneNucleonOvlap[];

extern ManyBodyOperators OpOneNucleonOvlaps;



void initOpOneNucleonOvlaps(OneNucleonOvlapsPara* par);


void calcOneNucleonOvlapsod(OneNucleonOvlapsPara* par,
			    const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    complex double* specamplitude);

void writeOneNucleonOvlaps(FILE* fp,
			   const Projection* P,
			   const OneNucleonOvlapsPara* p,
			   int jfin, int pfin, int afin,
			   int jini, int pini, int aini,
			   void* specamplitudeme,
			   const Eigenstates* Efin,
			   const Eigenstates* Eini);

#endif

