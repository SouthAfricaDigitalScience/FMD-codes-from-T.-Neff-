/**

  \file RadiiAll.h

  calculate radii


  (c) 2009 Thomas Neff

*/


#ifndef _RADIIALL_H
#define _RADIIALL_H

#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"

typedef struct {
  double r2m;
  double r2p;
  double r2n;
  double r2pp;
  double r2nn;
} RadiiAll;

typedef struct {
  complex double r2m;
  complex double r2p;
  complex double r2n;
  complex double r2pp;
  complex double r2nn;
} RadiiAllod;


extern ManyBodyOperator OpRadiiAll;


void calcRadii2all(const SlaterDet* Q, const SlaterDetAux* X,
		   RadiiAll* r2);


void calcRadii2allod(void* par,
		     const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X,
		     RadiiAllod* r2);

void writeprojectedRadiiAll(FILE* fp,
			    const SlaterDet* Q,
			    const Projection* P,
			    RadiiAllod** radiiall,
			    const Eigenstates* E);

#endif
