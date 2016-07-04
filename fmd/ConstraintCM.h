/**

  \file ConstraintCM.h

  constrain CM position in coordinate and momentum space


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINCM_H
#define _CONSTRAINCM_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintCM;


void calcConstraintCM(const SlaterDet* Q, const SlaterDetAux* X,
		      double* cmconstraint);

void calcgradConstraintCM(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

double outputConstraintCM(double val);

#endif
