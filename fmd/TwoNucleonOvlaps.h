/**

   \file TwoNucleonOvlaps.h

   calculate two-nucleon spectroscopic amplitudes in momentum space


   (c) 2012 Thomas Neff

*/


#ifndef _TWONUCLEONOVLAPS_H
#define _TWONUCLEONOVLAPS_H


#include <stdio.h>

#include "SlaterDet.h"
#include "Projection.h"


typedef struct {
  double qmax;
  int npoints;
  int nalpha;
  int ncosb;
  int recoil;
} TwoNucleonOvlapsPara;

#define NSPECT 5
#define DIMSPECT 5

#define NSPECY 5
#define DIMSPECY 5

extern ManyBodyOperator OpTwoNucleonOvlapT[];
extern ManyBodyOperator OpTwoNucleonOvlapY[];

extern ManyBodyOperators OpTwoNucleonOvlapsT;
extern ManyBodyOperators OpTwoNucleonOvlapsY;


void initOpTwoNucleonOvlapsT(TwoNucleonOvlapsPara* par);

void initOpTwoNucleonOvlapsY(TwoNucleonOvlapsPara* par);


void calcTwoNucleonOvlapsTod(TwoNucleonOvlapsPara* par,
			     const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     complex double* specamplitude);

void calcTwoNucleonOvlapsYod(TwoNucleonOvlapsPara* par,
			     const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     complex double* specamplitude);


void writeTwoNucleonOvlapsT(FILE* fp,
			    const Projection* P,
			    const TwoNucleonOvlapsPara* p,
			    int jfin, int pfin, int afin,
			    int jini, int pini, int aini,
			    void* specamplitudeme,
			    const Eigenstates* Efin,
			    const Eigenstates* Eini);

void writeTwoNucleonOvlapsY(FILE* fp,
			    const Projection* P,
			    const TwoNucleonOvlapsPara* p,
			    int jfin, int pfin, int afin,
			    int jini, int pini, int aini,
			    void* specamplitudeme,
			    const Eigenstates* Efin,
			    const Eigenstates* Eini);

#endif

