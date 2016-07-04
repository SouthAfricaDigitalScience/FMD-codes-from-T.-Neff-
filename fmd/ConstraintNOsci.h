/**

  \file ConstraintNOsci.h

  Constrain number of oscillator quanta
  for an isotropic oscillator


  (c) 2006 Thomas Neff

*/


#ifndef _CONSTRAINTNOSCI_H
#define _CONSTRAINTNOSCI_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"

extern Constraint ConstraintNOsci;
extern Constraint ConstraintPNOsci;
extern Constraint ConstraintNNOsci;


void calcConstraintNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
			 double* nosci);

void calcConstraintPNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
                          double* nosci);

void calcConstraintNNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
                          double* nosci);


void calcgradConstraintNOsci(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     gradSlaterDet* dnosci);

void calcgradConstraintPNOsci(const SlaterDet* Q, const SlaterDetAux* X,
                              const gradSlaterDetAux* dX,
                              gradSlaterDet* dnosci);

void calcgradConstraintNNOsci(const SlaterDet* Q, const SlaterDetAux* X,
                              const gradSlaterDetAux* dX,
                              gradSlaterDet* dnosci);


double outputConstraintNOsci(double val);


#endif


