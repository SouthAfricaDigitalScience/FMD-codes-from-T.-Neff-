/**

  \file Potential.h

  calculate matrix element of Interaction


  (c) 2003 Thomas Neff

*/


#ifndef _POTENTIAL_H
#define _POTENTIAL_H


#include <complex.h>

#include "SlaterDet.h"
#include "Interaction.h"


void calcPotential(const Interaction *Int,
		   const SlaterDet* Q, const SlaterDetAux* X, double v[]);

void calcPotentialrowcol(const Interaction *Int,
			 const SlaterDet* Q, const SlaterDetAux* X, double v[],
			 int k, int l);

void calcPotentialod(const Interaction *Int,
		     const SlaterDet* Q, const SlaterDet* Qp, 
		     const SlaterDetAux* X, complex double v[]);

void calcPotentialodrowcol(const Interaction *Int,
			   const SlaterDet* Q, const SlaterDet* Qp, 
			   const SlaterDetAux* X, complex double v[],
			   int k, int l);

void calcPotentialsHF(const Interaction* Int,
		      const SlaterDet* Q, const SlaterDetAux* X,
		      void* mes);

void calcPotentialHF(const Interaction* Int,
		     const SlaterDet* Q, const SlaterDetAux* X,
		     complex double* mes);

#endif
