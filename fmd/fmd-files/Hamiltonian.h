/**

  \file Hamiltonian.h

  calculate matrix element of Hamiltonian


  (c) 2003 Thomas Neff

*/


#ifndef _HAMILTONIAN_H
#define _HAMILTONIAN_H


#include <complex.h>

#include "SlaterDet.h"
#include "Interaction.h"


void calcHamiltonian(const Interaction *Int,
		     const SlaterDet* Q, const SlaterDetAux* X,
		     double* h);

void calcHamiltonianod(const Interaction *Int,
		       const SlaterDet* Q, const SlaterDet* Qp, 
		       const SlaterDetAux* X,
		       complex double* h);

void calcHamiltonianHF(const Interaction* Int,
		       const SlaterDet* Q, const SlaterDetAux* X,
		       void* mes);

#endif
