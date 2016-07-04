/**

  \file Oscillator.h

  external quadratic and quartic oscillator


  (c) 2003 Thomas Neff

*/


#ifndef _OSCILLATOR_H
#define _OSCILLATOR_H


#include "SlaterDet.h"


void calcOsci2(const SlaterDet* Q, const SlaterDetAux* X, 
	       double omega, double* hosci);


void calcOsci4(const SlaterDet* Q, const SlaterDetAux* X, 
	       double kappa, double* hosci);


#endif
