/**

  \file ConstraintQuadrupole.h

  Constrain mass quadrupole moment


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTQUADRUPOLE_H
#define _CONSTRAINTQUADRUPOLE_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"

extern Constraint ConstraintQ2;
extern Constraintod ConstraintQ2od;

extern Constraint ConstraintEQ2;
extern Constraintod ConstraintEQ2od;

extern Constraint ConstraintNQ2;
extern Constraintod ConstraintNQ2od;

extern Constraint ConstraintDetQ;

extern Constraint ConstraintQ2offdiagonal;


void calcConstraintQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* q2);

void calcgradConstraintQ2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dq2);

void calcConstraintQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double* q2);

void calcgradConstraintQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    const double* q,
			    gradSlaterDet* dq2);

double outputConstraintQ2(double val);


void calcConstraintEQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* q2);

void calcgradConstraintEQ2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dq2);

void calcConstraintEQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double* q2);

void calcgradConstraintEQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* q,
			     gradSlaterDet* dq2);

double outputConstraintEQ2(double val);


void calcConstraintNQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* q2);

void calcgradConstraintNQ2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dq2);

void calcConstraintNQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double* q2);

void calcgradConstraintNQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* q,
			     gradSlaterDet* dq2);

double outputConstraintNQ2(double val);


void calcConstraintDetQ(const SlaterDet* Q, const SlaterDetAux* X, 
			double* detq);


void calcgradConstraintDetQ(const SlaterDet* Q, const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* ddetq);

double outputConstraintDetQ(double val);


void calcConstraintQ2offdiagonal(const SlaterDet* Q, const SlaterDetAux* X,
				 double* dq2);

void calcgradConstraintQ2offdiagonal(const SlaterDet* Q, const SlaterDetAux* X,
				     const gradSlaterDetAux* dX,
				     gradSlaterDet* dq2);

double outputConstraintQ2offdiagonal(double val);

#endif
