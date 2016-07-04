/**

  \file TimeReversal.h

  calculate TimeReversal


  (c) 2007 Thomas Neff

*/


#ifndef _TIMEREVERSAL_H
#define _TIMEREVERSAL_H

#include "SlaterDet.h"
#include "Projection.h"


extern ManyBodyOperator OpTimeReversal;


void calcTimeReversal(const SlaterDet* Q, const SlaterDetAux* X, 
		      complex double* t);

void calcTimeReversalod(void* dummy,
			const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X, 
			complex double* t);

void showprojectedTimeReversal(FILE* fp,
			       const Projection* P,
			       complex double** tr,
			       const Eigenstates* E);


#endif
