/**

  \file ParameterizationFMD.h

  Parametrization of SlaterDet.

  One-to-one mapping between Parameterization and SlaterDet


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMD_H
#define _PARAMETERIZATIONFMD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMD;


int FMDread(FILE* fp, Para* q);
int FMDwrite(FILE* fp, const Para* q);

void FMDclone(const Para* q, Para* qp);

void FMDinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDtoSlaterDet(const Para* q, SlaterDet* Q);
void FMDprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 

void SlaterDetinitFMD(const SlaterDet* Q, Para* q);

#endif
