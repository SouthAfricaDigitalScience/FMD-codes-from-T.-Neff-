/**

  \file NOscillator.h

  single-particle harmonic oscillator hamiltonian


  (c) 2005 Thomas Neff

*/


#ifndef _NOSCILLATOR_H
#define _NOSCILLATOR_H

#include "SlaterDet.h"


void calcHOsciHF(const SlaterDet* Q, const SlaterDetAux* X, double omega,
		 void* mes);

#endif
