/**

  \file gradOscillator.h

  calculate gradient of kinetic energy


  (c) 2003 Thomas Neff

*/


#ifndef _GRADOSCILLATOR_H
#define _GRADOSCILLATOR_H

#include "SlaterDet.h"
#include "gradSlaterDet.h"


void calcgradOsci2(const SlaterDet* Q, const SlaterDetAux* X, 
		   const gradSlaterDetAux* dX,
		   double omega, gradSlaterDet* dhosci);


void calcgradOsci4(const SlaterDet* Q, const SlaterDetAux* X, 
		   const gradSlaterDetAux* dX,
		   double kappa, gradSlaterDet* dhosci);


#endif
