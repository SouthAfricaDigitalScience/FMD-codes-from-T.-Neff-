/**

  \file Parity.h

  calculate parity


  (c) 2003 Thomas Neff

*/


#ifndef _PARITY_H
#define _PARITY_H

#include "SlaterDet.h"


void calcParity(const SlaterDet* Q, const SlaterDetAux* X, 
		double* pi);

void calcParityod(const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double* pi);

void calcparHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes);

#endif
