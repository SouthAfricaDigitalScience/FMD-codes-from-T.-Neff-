/**

  \file Constraintl2.c

  Constrain expectation value of l2


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTl2_H
#define _CONSTRAINTl2_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint Constraintl2;


void calcConstraintl2(const SlaterDet* Q, const SlaterDetAux* X,
		      double* l2);


void calcgradConstraintl2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dl2);

double outputConstraintl2(double val);

#endif
