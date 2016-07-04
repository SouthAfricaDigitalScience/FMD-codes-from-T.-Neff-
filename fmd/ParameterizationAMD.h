/**

  \file ParameterizationAMD.h

  Parametrization of SlaterDet.

  common real width for all single-particle states


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONAMD_H
#define _PARAMETERIZATIONAMD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAMD;


int AMDread(FILE* fp, Para* q);
int AMDwrite(FILE* fp, const Para* q);

void AMDclone(const Para* q, Para* qp);

void AMDinitSlaterDet(const Para* q, SlaterDet* Q);
void AMDtoSlaterDet(const Para* q, SlaterDet* Q);
void AMDprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 


#endif
