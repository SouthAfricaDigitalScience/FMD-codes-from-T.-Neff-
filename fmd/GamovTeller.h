/**

  \file GamovTeller.h 

  GamovTeller transitions

  calculates Gamov Teller transition matrix elements
  in spherical basis


  (c) 2004 Thomas Neff

*/


#ifndef _GAMOVTELLER_H
#define _GAMOVTELLER_H


#include "SlaterDet.h"
#include "Projection.h"


// GTplus sigma tau+
extern ManyBodyOperator OpGTplus;

// GTminus sigma tau-
extern ManyBodyOperator OpGTminus;



void calcGTplusod(void* Par,
		  const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double GTplus[3]);


void calcGTminusod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double GTminus[3]);

					 
void writeprojectedtransitionGTplus(FILE* fp,
				    const Projection* P,
				    const complex double**** gtp,
				    const Eigenstates* E, 
				    const Eigenstates* Ep);


#endif
