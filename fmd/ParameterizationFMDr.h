/**

  \file ParameterizationFMDr.h

  Parametrization of SlaterDet.

  width parameter a purely real


  (c) 2007 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMDR_H
#define _PARAMETERIZATIONFMDR_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMDr;


int FMDrread(FILE* fp, Para* q);
int FMDrwrite(FILE* fp, const Para* q);

void FMDrclone(const Para* q, Para* qp);

void FMDrinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDrtoSlaterDet(const Para* q, SlaterDet* Q);
void FMDrprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 

void SlaterDetinitFMDr(const SlaterDet* Q, Para* q);

#endif
