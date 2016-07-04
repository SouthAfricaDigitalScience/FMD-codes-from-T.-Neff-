/**

  \file gradHamiltonian.h

  calculate gradients of Hamiltonian matrix element


  (c) 2003 Thomas Neff

*/


#ifndef _GRADHAMILTONIAN_H
#define _GRADHAMILTONIAN_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Interaction.h"


/// calculates gradient of Hamiltonian
/// overwrites gradient dh
void calcgradHamiltonian(const Interaction *P,
			 const SlaterDet* Q, const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dh);

void calcgradHamiltonianod(const Interaction *P,
			   const SlaterDet* Q, const SlaterDet* Qp, 
			   const SlaterDetAux* X, 
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dh);


#endif
