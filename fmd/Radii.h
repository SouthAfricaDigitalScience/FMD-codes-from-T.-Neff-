/**

  \file Radii.h

  calculate radii


  (c) 2003 Thomas Neff

*/


#ifndef _RADII_H
#define _RADII_H

#include "SlaterDet.h"


void calcRadii2(const SlaterDet* Q, const SlaterDetAux* X, 
		double* r2mass, double* r2proton, double* r2neutron);

void calcRadii2od(const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double* r2mass,
		  complex double* r2proton, complex double* r2neutron);

#endif
