/**

  \file ParameterizationAlphaA.h

  Parametrization of SlaterDet.

  alpha-cluster nuclei with fixed width a 


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONALPHAA_H
#define _PARAMETERIZATIONALPHAA_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAlphaA;

int AlphaAread(FILE* fp, Para* q);
int AlphaAwrite(FILE* fp, const Para* q);

void AlphaAclone(const Para* q, Para* qp);

void AlphaAinitSlaterDet(const Para* q, SlaterDet* Q);
void AlphaAtoSlaterDet(const Para* q, SlaterDet* Q);
void AlphaAprojectgradSlaterDet(const Para* q, 
			     const gradSlaterDet* dQ, double* dq); 


#endif
