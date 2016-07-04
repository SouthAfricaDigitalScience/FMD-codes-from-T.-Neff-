/**

  \file gradKineticEnergy.h

  calculate gradient of kinetic energy


  (c) 2003 Thomas Neff

*/


#ifndef _GRADKINETICENERGY_H
#define _GRADKINETICENERGY_H

#include "SlaterDet.h"
#include "gradSlaterDet.h"

void calcgradT(const SlaterDet* Q, const SlaterDetAux* X, 
	       const gradSlaterDetAux* dX,
	       gradSlaterDet* dt);

void calcgradTod(const SlaterDet* Q, const SlaterDet* Qp, 
		 const SlaterDetAux* X, 
		 const gradSlaterDetAux* dX,
		 gradSlaterDet* dt);


#endif
