/**

  \file CenterofMass.h

  calculate observables related to center of mass motion


  (c) 2003 Thomas Neff

*/


#ifndef _CENTEROFMASS_H
#define _CENTEROFMASS_H

#include "SlaterDet.h"


void calcCMPosition(const SlaterDet* Q, const SlaterDetAux* X, 
		    double xcm[3]);

void calcCMVelocity(const SlaterDet* Q, const SlaterDetAux* X, 
		    double xcm[3]);


void calcTCM(const SlaterDet* Q, const SlaterDetAux* X, 
	     double* tcm);

void calcTCMod(const SlaterDet* Q, const SlaterDet* Qp, 
	       const SlaterDetAux* X, 
	       complex double* tcm);

void calcTCMHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes);


#endif
