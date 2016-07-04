/**

  \file SpatialOrientation.h

  determinas euler angles of spatial orientation
  by diagonalizing the inertia tensor


  (c) 2003 Thomas Neff

*/


#ifndef _SPATIALORIENTATION_H
#define _SPATIALORIENTATION_H

#include "SlaterDet.h"


/// calculate mass inertia tensor 
void calcInertiaTensor(const SlaterDet* Q, const SlaterDetAux* X,
		       double T[9]);


/// calculate euler angles to describe orientation of mass distribution
/// Slater determinant center of mass should be in origin
void calcSpatialOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			    double* alpha, double* beta, double* gamma);


void calcOrientedOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			     double* alpha, double* beta, double* gamma);


void calcSpinOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			 double* alpha, double* beta, double* gamma);


#endif
