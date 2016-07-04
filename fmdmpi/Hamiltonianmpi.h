/**

  \file Hamiltonianmpi.h

  calculate matrix element of Interaction


  (c) 2003 Thomas Neff

*/


#ifndef _HAMILTONIANMPI_H
#define _HAMILTONIANMPI_H


#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Interaction.h"


void calcHamiltonianmpi(const Interaction *P,
			const SlaterDet* Q, const SlaterDetAux* X,
			double* h);

void calcHamiltonianodmpi(const Interaction *P,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X,
			  complex double* h);

void calcPotentialmpi(const Interaction *P,
		      const SlaterDet* Q, const SlaterDetAux* X, double v[]);

void calcPotentialodmpi(const Interaction *P,
			const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X, complex double v[]);

#endif
