/**

  \file ParameterizationFMDS.h

  Parametrization of SlaterDet.

  FMD Parameterization of Slater determinant with fixed spins


  (c) 2004 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMDS_H
#define _PARAMETERIZATIONFMDS_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMDS;


int FMDSread(FILE* fp, Para* q);
int FMDSwrite(FILE* fp, const Para* q);

void FMDSclone(const Para* q, Para* qp);

void FMDSinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDStoSlaterDet(const Para* q, SlaterDet* Q);
void FMDSprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 

#endif
