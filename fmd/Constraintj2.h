/**

  \file ConstraintJ2.c

  Constrain expectation value of J2


  (c) 2003 Thomas Neff

*/


#ifndef _ConstraintJ2_H
#define _ConstraintJ2_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintJ2;

extern Constraintod ConstraintJ2od;


void calcConstraintJ2(const SlaterDet* Q, const SlaterDetAux* X,
		      double* j2);

void calcgradConstraintJ2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dj2);

void calcConstraintJ2od(const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X,
			complex double* j2);

void calcgradConstraintJ2od(const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dj2);

double outputConstraintJ2(double val);

#endif
