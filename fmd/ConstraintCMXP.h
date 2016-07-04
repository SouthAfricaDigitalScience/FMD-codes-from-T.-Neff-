/**

  \file ConstraintCM.h

  constrain CM position in coordinate and momentum space


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTCMXP_H
#define _CONSTRAINTCMXP_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintX;
extern Constraint ConstraintY;
extern Constraint ConstraintZ;

extern Constraint ConstraintPX;
extern Constraint ConstraintPY;
extern Constraint ConstraintPZ;


void calcConstraintX(const SlaterDet* Q, const SlaterDetAux* X,
		     double* cmconstraint);

void calcConstraintY(const SlaterDet* Q, const SlaterDetAux* X,
		     double* cmconstraint);

void calcConstraintZ(const SlaterDet* Q, const SlaterDetAux* X,
		     double* cmconstraint);

void calcConstraintPX(const SlaterDet* Q, const SlaterDetAux* X,
		      double* cmconstraint);

void calcConstraintPY(const SlaterDet* Q, const SlaterDetAux* X,
		      double* cmconstraint);

void calcConstraintPZ(const SlaterDet* Q, const SlaterDetAux* X,
		      double* cmconstraint);

void calcgradConstraintX(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

void calcgradConstraintY(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

void calcgradConstraintZ(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

void calcgradConstraintPX(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

void calcgradConstraintPY(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

void calcgradConstraintPZ(const SlaterDet* Q, const SlaterDetAux* X,
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dcmconstraint);

double outputConstraintCM(double val);

#endif
