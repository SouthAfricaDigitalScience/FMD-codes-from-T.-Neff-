/**

   \file MeanOscillatorQuanta.h

   calculate mean number of oscillator quanta in many-body state


   (c) 2007 Thomas Neff

*/


#ifndef _MEANOSCILLATORQUANTA_H
#define _MEANOSCILLATORQUANTA_H

#include "SlaterDet.h"
#include "Projection.h"


typedef struct {
  complex double r2[3];
  complex double p2[3];
} MeanOsciQuanta;
 

// r2 and p2 operators
extern ManyBodyOperator OpMeanOsciQuanta;

void calcMeanOsciQuantaod(void* Par,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X,
			  MeanOsciQuanta* q);

void writeprojectedMeanOscillatorQuanta(FILE* fp, double omega,
					const Projection* P,
					const MeanOsciQuanta** Q,
					const Eigenstates* E);

#endif
