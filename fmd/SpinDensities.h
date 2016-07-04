/**

  \file SpinDensities.h

  calculate one-body spin densities in coordinate and momentum space


  (c) 2003 Thomas Neff

*/


#ifndef _SPINDENSITIES_H
#define _SPINDENSITIES_H


void calcSpinDensitiesCoordinate(const SlaterDet* Q, const SlaterDetAux* X, 
				 int v, int n, double x, double* dens);

void calcSpinDensitiesMomentum(const SlaterDet* Q, const SlaterDetAux* X, 
			       int v, int n, double p, double* dens);


#endif


