/**

  \file ConstraintOctupole.h

  Constrain Mass Octupole
  center of mass corrected by substraction of Xcm expectation value


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINTOCTUPOLE_H
#define _CONSTRAINTOCTUPOLE_H

#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Constraint.h"


extern Constraint ConstraintO2;

extern Constraint ConstraintEO2;

extern Constraint ConstraintNO2;


void calcConstraintO2(const SlaterDet* Q, const SlaterDetAux* X,
		      double* o2);

void calcConstraintEO2(const SlaterDet* Q, const SlaterDetAux* X,
		       double* o2);

void calcConstraintNO2(const SlaterDet* Q, const SlaterDetAux* X,
		       double* o2);

void calcgradConstraintO2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* do2);

void calcgradConstraintEO2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* do2);

void calcgradConstraintNO2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* do2);

double outputConstraintO2(double val);

double outputConstraintEO2(double val);

double outputConstraintNO2(double val);


#endif


