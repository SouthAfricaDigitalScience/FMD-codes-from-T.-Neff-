/**

  \file Observables.h

  calculate observables for Slater determinants


  (c) 2003 Thomas Neff

*/


#ifndef _OBSERVABLES_H
#define _OBSERVABLES_H

#include <stdio.h>

#include "SlaterDet.h"
#include "Interaction.h"


typedef struct {
  double h;
  double tcm;
  double r2m;
  double r2p;
  double r2n;
  double l2;
  double s2;
  double j2;
  double pi;
  double t2;
  double t;
  double v[MAXINTERACTIONS+1];
} Observables;


typedef struct {
  complex double n;
  complex double h;
  complex double tcm;
  complex double r2m;
  complex double r2p;
  complex double r2n;
  complex double l2;
  complex double s2;
  complex double j2;
  complex double pi;
  complex double t2;
  complex double t;
  complex double v[MAXINTERACTIONS+1];
} Observablesod;



void calcObservables(const Interaction* Int,
		     const SlaterDet* Q,
		     Observables* obs);


void calcObservablesParity(const Interaction* Int, int parity,
			   const SlaterDet* Q, Observables* obs);


void fprintObservables(FILE* fp, 
		       const Interaction* Int, const SlaterDet* Q,
		       const Observables* obs);


void calcObservablesod(const Interaction* Int,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       const SlaterDetAux* X,
		       Observablesod* obs);


void scaleObservablesod(const Interaction* Int, Observablesod* obs);


#endif
