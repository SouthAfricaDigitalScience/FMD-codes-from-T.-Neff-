/**

  \file Densities.h

  calculate one-body densities in coordinate and momentum space


  (c) 2003 Thomas Neff

*/


#ifndef _DENSITIES_H
#define _DENSITIES_H


void calcDensitiesCoordinate(const SlaterDet* Q,
			     int v, int n, double x, double* dens);

void calcDensitiesCoordinateParity(const SlaterDet* Q, int par,
				   int v, int n, double x, double* dens);


void calcRadialDensitiesCoordinate(const SlaterDet* Q,
				   int izcw, int n, double r, double* dens);


void calcDensitiesMomentum(const SlaterDet* Q,
			   int v, int n, double p, double* dens);


void calcRadialDensitiesMomentum(const SlaterDet* Q,
				 int izcw, int n, double p, double* dens);


void calcDensitiesCoordinateHF(const SlaterDet* Q, const SlaterDetAux* X, 
                               int v, int n, double x, void* mes);


void calcDensitiesMomentumHF(const SlaterDet* Q, const SlaterDetAux* X, 
                             int v, int n, double x, void* mes);


#endif


