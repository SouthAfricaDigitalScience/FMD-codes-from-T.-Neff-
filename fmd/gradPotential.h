/**

  \file gradPotenial.h

  calculate gradients of interaction matrix element


  (c) 2003 Thomas Neff

*/


#ifndef _GRADPOTENTIAL_H
#define _GRADPOTENTIAL_H


#include <complex.h>

#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Interaction.h"


void calcgradPotential(const Interaction *P,
		       const SlaterDet* Q, const SlaterDetAux* X, 
		       const gradSlaterDetAux* dX,
		       gradSlaterDet* dv);

void calcgradPotentialod(const Interaction *P,
			 const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dv);

void calcgradPotentialrowcol(const Interaction *P,
			     const SlaterDet* Q, const SlaterDetAux* X, 
			     const gradSlaterDetAux* dX,
			     gradSlaterDet* dv,
			     int k, int l);

void calcgradPotentialodrowcol(const Interaction *P,
			       const SlaterDet* Q, const SlaterDet* Qp,
			       const SlaterDetAux* X, 
			       const gradSlaterDetAux* dX,
			       gradSlaterDet* dv,
			       int k, int l);

#endif
