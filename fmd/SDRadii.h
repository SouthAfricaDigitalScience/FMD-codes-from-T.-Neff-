/**

  \file SDRadii.h

  calculate radii without correcting for center of mass


  (c) 2003 Thomas Neff

*/


#ifndef _SDRADII_H
#define _SDRADII_H

#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"

typedef struct {
  double r2m;
  double r2p;
  double r2n;
} SDRadii;

typedef struct {
  complex double r2m;
  complex double r2p;
  complex double r2n;
} SDRadiiod;


extern ManyBodyOperator OpSDRadii;


void calcSDRadii2(const SlaterDet* Q, const SlaterDetAux* X,
		  SDRadii* r2);

void calcSDRadii2od(void* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, 
		    SDRadiiod* r2);

void writeprojectedSDRadii(FILE* fp,
			   const SlaterDet* Q,
			   const Projection* P,
			   SDRadiiod** radii,
			   const Eigenstates* E);

#endif
