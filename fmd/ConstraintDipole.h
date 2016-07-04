/**

  \file ConstraintDipole.h

  Constrain dipole moment


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTEDIPOLE_H
#define _CONSTRAINTEDPOLE_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"

extern Constraint ConstraintED2;


void calcConstraintED2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* ed2);


void calcgradConstraintED2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* ded2);

double outputConstraintED2(double val);


#endif
