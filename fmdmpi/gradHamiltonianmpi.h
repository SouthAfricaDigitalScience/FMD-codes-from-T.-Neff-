/**

  \file gradHamiltonianmpi.h

  calculate gradients of Hamiltonian matrix element


  (c) 2003 Thomas Neff

*/


#ifndef _GRADHAMILTONIANMPI_H
#define _GRADHAMILTONIANMPI_H


#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Interaction.h"


/// calculates gradient of Hamiltonian
/// overwrites gradient dh
void calcgradHamiltonianmpi(const Interaction *P,
			    const SlaterDet* Q, const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* dh);

void calcgradHamiltonianodmpi(const Interaction *P,
			      const SlaterDet* Q, const SlaterDet* Qp,
			      const SlaterDetAux* X, 
			      const gradSlaterDetAux* dX,
			      gradSlaterDet* dh);

void calcgradPotentialmpi(const Interaction *P,
                          const SlaterDet* Q, const SlaterDetAux* X, 
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dv);

void calcgradPotentialodmpi(const Interaction *P,
			    const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* dv);

#endif
