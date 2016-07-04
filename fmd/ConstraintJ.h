/**

  \file ConstraintJ.c

  Constrain expectation values of Jx, Jy, Jz


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTJ_H
#define _CONSTRAINTJ_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"


void calcConstraintJ(const SlaterDet* Q, const SlaterDetAux* X,
		     double* j2);

void calcgradConstraintJ(const SlaterDet* Q, const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dj2);


#endif
