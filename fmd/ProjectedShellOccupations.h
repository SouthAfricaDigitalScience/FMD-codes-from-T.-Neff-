/**

   \file ProjectedShellOccupations.c

   calculate the many-body occupations in a harmonic oscillator basis
   see Suzuki et al, PRC 54, 2073 (1996)


   (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"


#ifndef _PROJECTEDSHELLOCCUPATIONS_H
#define _PROJECTEDSHELLOCCUPATIONS_H


typedef struct {
  double omega;         // oscillator parameter [fm^-1]
  int nmax;             // occupations up to nmax hw
  int nint;             // number of grid points for theta integration
} ShellOccupationsPara;


extern ManyBodyOperator OpShellOccupations;


void initOpShellOccupations(ShellOccupationsPara* par);

void calcShellOccupationsod(ShellOccupationsPara* par,
			    const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    complex double *shelloccupations);

void writeprojectedShellOccupations(FILE* fp, const ShellOccupationsPara* par,
				    const Projection* P,
				    void* shelloccupations,
				    const Eigenstates* E,
				    int totalonly);

#endif
