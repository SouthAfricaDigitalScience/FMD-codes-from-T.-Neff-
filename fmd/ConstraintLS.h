/**

  \file ConstraintLS.c

  Constrain expectation value of LS


  (c) 2005 Thomas Neff

*/


#ifndef _CONSTRAINTLS_H
#define _CONSTRAINTLS_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintLS;

extern Constraintod ConstraintLSod;


void calcConstraintLS(const SlaterDet* Q, const SlaterDetAux* X,
		      double* ls);

void calcgradConstraintLS(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dls);

void calcConstraintLSod(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X,
			complex double* ls);

void calcgradConstraintLSod(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dls);

double outputConstraintLS(double val);

#endif
