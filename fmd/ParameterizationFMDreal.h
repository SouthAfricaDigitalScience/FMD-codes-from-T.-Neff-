/**

  \file ParameterizationFMDreal.h

  Parametrization of SlaterDet.

  Parameterization only uses real parts of a and b and SlaterDet


  (c) 2005 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMDREAL_H
#define _PARAMETERIZATIONFMDREAL_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMDreal;


int FMDrealread(FILE* fp, Para* q);
int FMDrealwrite(FILE* fp, const Para* q);

void FMDrealclone(const Para* q, Para* qp);

void FMDrealinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDrealtoSlaterDet(const Para* q, SlaterDet* Q);
void FMDrealprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 

void SlaterDetinitFMDreal(const SlaterDet* Q, Para* q);

#endif
