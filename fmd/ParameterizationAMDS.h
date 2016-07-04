/**

  \file ParameterizationAMDS.h

  Parametrization of SlaterDet.

  common width for all single-particle states
  fixed spin


  (c) 2006 Thomas Neff

*/


#ifndef _PARAMETERIZATIONAMDS_H
#define _PARAMETERIZATIONAMDS_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAMDS;


int AMDSread(FILE* fp, Para* q);
int AMDSwrite(FILE* fp, const Para* q);

void AMDSclone(const Para* q, Para* qp);

void AMDSinitSlaterDet(const Para* q, SlaterDet* Q);
void AMDStoSlaterDet(const Para* q, SlaterDet* Q);
void AMDSprojectgradSlaterDet(const Para* q, 
			      const gradSlaterDet* dQ, double* dq); 


#endif
