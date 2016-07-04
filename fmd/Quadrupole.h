/**

  \file Quadrupole.c

  calculate quadrupole parameters beta and gamma
  for nucleons, protons and neutrons


  (c) 2003, 2005 Thomas Neff

*/


#ifndef _QUADRUPOLE_H
#define _QUADRUPOLE_H

#include "SlaterDet.h"


void calcQuadrupole(const SlaterDet* Q, const SlaterDetAux* X, 
		    double beta[3], double gamma[3]);

void calcQuadrupolePrime(const SlaterDet* Q, const SlaterDetAux* X, 
			 double beta[3]);


#endif
