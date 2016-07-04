/**

  \file ParameterizationFMDvxz.h

  Parametrization of SlaterDet with symmetry with respect to reflection on xz-plane


  (c) 2009 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMDVXZ_H
#define _PARAMETERIZATIONFMDVXZ_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMDvxz;


int FMDvxzread(FILE* fp, Para* q);
int FMDvxzwrite(FILE* fp, const Para* q);

void FMDvxzclone(const Para* q, Para* qp);

void FMDvxzinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDvxztoSlaterDet(const Para* q, SlaterDet* Q);
void FMDvxzprojectgradSlaterDet(const Para* q, 
                                const gradSlaterDet* dQ, double* dq); 

#endif
