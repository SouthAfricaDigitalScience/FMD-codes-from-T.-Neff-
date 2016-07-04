/**

  \file gradCenterofMass.h

  calculate gradient of center of mass energy


  (c) 2003 Thomas Neff

*/


#ifndef _GRADCENTEROFMASS_H
#define _GRADCENTEROFMASS_H

#include "SlaterDet.h"
#include "gradSlaterDet.h"


/// calculate the gradient of -Tcm !
void calcgradTCM(const SlaterDet* Q, const SlaterDetAux* X,
		 const gradSlaterDetAux* dX,
		 gradSlaterDet* G);

void calcgradTCMod(const SlaterDet* Q, const SlaterDet* Qp, 
		   const SlaterDetAux* X,
		   const gradSlaterDetAux* dX,
		   gradSlaterDet* G);

#endif
