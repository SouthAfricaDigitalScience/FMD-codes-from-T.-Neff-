/**

  \file ConstraintT2.h

  Constrain expectation value of T2


  (c) 2004 Thomas Neff

*/


#ifndef _CONSTRAINTT2_H
#define _CONSTRAINTT2_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintT2;
extern Constraint ConstraintPT2;
extern Constraint ConstraintNT2;

extern Constraintod ConstraintT2od;


void calcConstraintT2(const SlaterDet* Q, const SlaterDetAux* X,
		      double* t2);

void calcConstraintPT2(const SlaterDet* Q, const SlaterDetAux* X,
                       double* t2);

void calcConstraintNT2(const SlaterDet* Q, const SlaterDetAux* X,
                       double* t2);

void calcgradConstraintT2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dt2);

void calcgradConstraintPT2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dt2);

void calcgradConstraintNT2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dt2);


void calcConstraintT2od(const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X,
			complex double* t2);

void calcgradConstraintT2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dt2);

double outputConstraintT2(double val);

#endif
