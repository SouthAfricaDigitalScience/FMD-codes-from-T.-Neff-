/**

  \file Densities3d.h

  calculate one-body densities in coordinate space


  (c) 2003 Thomas Neff

*/


#ifndef _DENSITIES3D_H
#define _DENSITIES3D_H


void calcDensitiesCoordinate3d(const SlaterDet* Q, const SlaterDetAux* X, 
			       int n, double x, double* dens);

void calcDensitiesMomentum3d(const SlaterDet* Q, const SlaterDetAux* X, 
			     int n, double p, double* dens);


#endif


