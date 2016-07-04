/**

  \file ConstraintMainAxes.h

  Constrain to main axes of quadrupole tensor


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTMAINAXES_H
#define _CONSTRAINTMAINAXES_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"


void calcConstraintMainAxes(const SlaterDet* Q, const SlaterDetAux* X, 
			    double* q2);


void calcgradConstraintMainAxes(const SlaterDet* Q, const SlaterDetAux* X,
				const gradSlaterDetAux* dX,
				gradSlaterDet* dq2);


#endif
