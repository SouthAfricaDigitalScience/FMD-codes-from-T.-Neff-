/**

  \file ConstraintS2.h

  Constrain expectation value of S2


  (c) 2004 Thomas Neff

*/


#ifndef _CONSTRAINTS2_H
#define _CONSTRAINTS2_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintS2;
extern Constraint ConstraintPS2;
extern Constraint ConstraintNS2;

extern Constraintod ConstraintS2od;


void calcConstraintS2(const SlaterDet* Q, const SlaterDetAux* X,
		      double* s2);

void calcConstraintPS2(const SlaterDet* Q, const SlaterDetAux* X,
                       double* s2);

void calcConstraintNS2(const SlaterDet* Q, const SlaterDetAux* X,
                       double* s2);

void calcgradConstraintS2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* ds2);

void calcgradConstraintPS2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* ds2);

void calcgradConstraintNS2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* ds2);

void calcConstraintS2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X,
			complex double* s2);

void calcgradConstraintS2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* ds2);

double outputConstraintS2(double val);

#endif
