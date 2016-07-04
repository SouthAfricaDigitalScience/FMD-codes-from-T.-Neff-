/**

  \file KineticEnergy.h

  calculate kinetic energy


  (c) 2003 Thomas Neff

*/


#ifndef _KINETICENERGY_H
#define _KINETICENERGY_H

#include "SlaterDet.h"


void calcT(const SlaterDet* Q, const SlaterDetAux* X, 
	   double* t);

void calcTod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X, complex double* t);

void calcTHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes);


#endif
