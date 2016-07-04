/**

  \file NOsci.h 

  Harmonic Oscillator quanta

  (c) 2007 Thomas Neff

*/


#ifndef _NOSCI_H
#define _NOSCI_H


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"

#include "numerics/cmath.h"
#include "misc/utils.h"
#include "misc/physics.h"


// NOsci
extern ManyBodyOperator OpNOsci;


void calcNOsciod(void* Par,
                 const SlaterDet* Q, const SlaterDet* Qp,
                 const SlaterDetAux* X, 
                 complex double nosci[2]);

void writeprojectedNOsci(FILE* fp,
                         const Projection* P,
                         const complex double (**nosci)[2],
                         const Eigenstates* E);


#endif
