/**

  \file ParameterizationAMDAS.h

  Parametrization of SlaterDet.

  common fixed width for all single-particle states
  fixed spin


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONAMDAS_H
#define _PARAMETERIZATIONAMDAS_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAMDAS;


int AMDASread(FILE* fp, Para* q);
int AMDASwrite(FILE* fp, const Para* q);

void AMDASclone(const Para* q, Para* qp);

void AMDASinitSlaterDet(const Para* q, SlaterDet* Q);
void AMDAStoSlaterDet(const Para* q, SlaterDet* Q);
void AMDASprojectgradSlaterDet(const Para* q, 
		    	       const gradSlaterDet* dQ, double* dq); 


#endif
