/**

  \file AngularMomenta.h

  calculate angular momenta


  (c) 2003 Thomas Neff

*/


#ifndef _ANGULARMOMENTA_H
#define _ANGULARMOMENTA_H

#include "SlaterDet.h"


/// calculate l2, s2 and j2
void calcAngularMomenta(const SlaterDet* Q, const SlaterDetAux* X, 
			double* l2, double* s2, double* j2);

void calcAngularMomentaod(const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X, 
			  complex double* l2, complex double* s2, 
			  complex double* j2);

/// calc single-particle HF angular momenta
void calcl2HF(const SlaterDet* Q, const SlaterDetAux* X,
	      void* mes);

void calcj2HF(const SlaterDet* Q, const SlaterDetAux* X,
	      void* mes);


/// calculate Jx^2, Jy^2, Jz^2
void calcJ2(const SlaterDet* Q, const SlaterDetAux* X, 
	    double j2[3]);

void calcJ2od(const SlaterDet* Q, const SlaterDet* Qp,
	      const SlaterDetAux* X, 
	      complex double j2[3]);


/// calculate Lx, Ly, Lz

void calcL(const SlaterDet* Q, const SlaterDetAux* X, 
	   double l[3]);

void calcLod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X, 
	     complex double l[3]);


/// calculate Sx, Sy, Sz

void calcS(const SlaterDet* Q, const SlaterDetAux* X, 
	   double s[3]);

void calcSod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X, 
	     complex double s[3]);

/// calculate Jx, Jy, Jz
void calcJ(const SlaterDet* Q, const SlaterDetAux* X, 
	   double j[3]);

void calcJod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X, 
	     complex double j[3]);

#endif
