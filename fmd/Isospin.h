/**

  \file Isospin.h

  calculate T^2


  (c) 2003 Thomas Neff

*/


#ifndef _ISOSPIN_H
#define _ISOSPIN_H

#include "SlaterDet.h"


void calcIsospin(const SlaterDet* Q, const SlaterDetAux* X, 
		 double* t2);

void calcIsospinod(const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double* t2);

void calct3HF(const SlaterDet* Q, const SlaterDetAux* X, void* mes);


#endif
