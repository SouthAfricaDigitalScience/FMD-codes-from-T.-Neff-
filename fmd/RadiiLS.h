/**

  \file RadiiLS.h

  calculate spin-orbit contribution to charge radii


  (c) 2011 Thomas Neff

*/


#ifndef _RADIILS_H
#define _RADIILS_H

#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"

typedef struct {
  double r2lsp;
  double r2lsn;
} RadiiLS;

typedef struct {
  complex double r2lsp;
  complex double r2lsn;
} RadiiLSod;


extern ManyBodyOperator OpRadiiLS;


void calcRadii2LS(const SlaterDet* Q, const SlaterDetAux* X,
		  RadiiLS* r2ls);


void calcRadii2LSod(void* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X,
		    RadiiLSod* r2ls);

void writeprojectedRadiiLS(FILE* fp,
			   const SlaterDet* Q,
			   const Projection* P,
			   const RadiiLSod** radiils,
			   const Eigenstates* E);

#endif
