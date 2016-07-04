/**

  \file ParameterizationAlphaC.h

  Parametrization of SlaterDet.

  alpha-cluster nuclei with fixed width a and vanishing momenta


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONALPHAC_H
#define _PARAMETERIZATIONALPHAC_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAlphaC;

int AlphaCread(FILE* fp, Para* q);
int AlphaCwrite(FILE* fp, const Para* q);

void AlphaCclone(const Para* q, Para* qp);

void AlphaCinitSlaterDet(const Para* q, SlaterDet* Q);
void AlphaCtoSlaterDet(const Para* q, SlaterDet* Q);
void AlphaCprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 


#endif
