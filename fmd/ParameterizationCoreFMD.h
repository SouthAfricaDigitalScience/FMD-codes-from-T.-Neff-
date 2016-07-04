/**

  \file ParameterizationCoreFMD.h

  Parametrization of SlaterDet.

  Core and FMD nucleons


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONCOREFMD_H
#define _PARAMETERIZATIONCOREFMD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationCoreFMD;


int CoreFMDread(FILE* fp, Para* q);
int CoreFMDwrite(FILE* fp, const Para* q);

void CoreFMDclone(const Para* q, Para* qp);

void CoreFMDinitSlaterDet(const Para* q, SlaterDet* Q);
void CoreFMDtoSlaterDet(const Para* q, SlaterDet* Q);
void CoreFMDprojectgradSlaterDet(const Para* q, 
				 const gradSlaterDet* dQ, double* dq); 

#endif
