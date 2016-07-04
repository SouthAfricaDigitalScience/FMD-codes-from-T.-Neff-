/**

  \file ElectroMagneticMultipole.c 

  Electro-Magnetic Multipole Operators

  calculates electric monopole, magnetic dipole and
  electric quadrupole moment in spherical basis


  (c) 2003 Thomas Neff

*/


#ifndef _EMMULTIPOLE_H
#define _EMMULTIPOLE_H

#include "SlaterDet.h"
#include "Projection.h"


/// Electric Monopole operator
extern ManyBodyOperator OpEMonopole;

/// Electric Dipole operator
extern ManyBodyOperator OpEDipole;

/// Magnetic Dipole operator
extern ManyBodyOperator OpMDipole;

/// Electric Quadrupole operator
extern ManyBodyOperator OpEQuadrupole;



void calcEMonopoleod(void* Par,
		     const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, 
		     complex double* emonopole);


void calcEDipoleod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double edipole[3]);


void calcMDipoleod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double mdipole[3]);


void calcEQuadrupoleod(void* Par,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       const SlaterDetAux* X, 
		       complex double equadrupole[5]);


void writeprojectedEMMultipoles(FILE* fp,
				const Projection* P,
				const complex double** emonopole,
				const complex double** edipole,
				const complex double** mdipole,
				const complex double** equadrupole,
				const Eigenstates* E);


void writeprojectedtransitionEMMultipoles(FILE* fp,
					  const Projection* P,
					  const complex double**** emonopole,
					  const complex double**** edipole,
					  const complex double**** mdipole,
					  const complex double**** equadrupole,
					  const Eigenstates* E);

#endif
