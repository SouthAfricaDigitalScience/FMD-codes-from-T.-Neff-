/**

  \file ParameterizationAMDd.h

  Parametrization of SlaterDet.

  two common real width parameters for all single-particle states


  (c) 2007 Thomas Neff

*/


#ifndef _PARAMETERIZATIONAMDD_H
#define _PARAMETERIZATIONAMDD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAMDd;


int AMDdread(FILE* fp, Para* q);
int AMDdwrite(FILE* fp, const Para* q);

void AMDdclone(const Para* q, Para* qp);

void AMDdinitSlaterDet(const Para* q, SlaterDet* Q);
void AMDdtoSlaterDet(const Para* q, SlaterDet* Q);
void AMDdprojectgradSlaterDet(const Para* q, 
			      const gradSlaterDet* dQ, double* dq); 


#endif
