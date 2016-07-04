/**

  \file ConstraintR2.h

  Constrain Mass and Charge Radius


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTR2_H
#define _CONSTRAINTR2_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"

extern Constraint ConstraintR2;

extern Constraintod ConstraintR2od;

extern Constraint ConstraintER2;

extern Constraintod ConstraintER2od;

extern Constraint ConstraintNR2;

extern Constraintod ConstraintNR2od;


void calcConstraintR2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* r2);


void calcgradConstraintR2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dr2);


void calcConstraintR2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double* r2);


void calcgradConstraintR2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dr2);


double outputConstraintR2(double val);


void calcConstraintER2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* r2);


void calcgradConstraintER2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dr2);


void calcConstraintER2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double* r2);


void calcgradConstraintER2od(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* dummy,
			     gradSlaterDet* dr2);


double outputConstraintER2(double val);


void calcConstraintNR2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* r2);


void calcgradConstraintNR2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dr2);


void calcConstraintNR2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double* r2);


void calcgradConstraintNR2od(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* dummy,
			     gradSlaterDet* dr2);


double outputConstraintNR2(double val);


#endif
