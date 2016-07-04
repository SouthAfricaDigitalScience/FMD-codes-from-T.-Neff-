/**

  \file ParameterizationAlpha.h

  Parametrization of SlaterDet.

  alpha-cluster nuclei


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONALPHA_H
#define _PARAMETERIZATIONALPHA_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAlpha;


int Alpharead(FILE* fp, Para* q);
int Alphawrite(FILE* fp, const Para* q);

void Alphaclone(const Para* q, Para* qp);

void AlphainitSlaterDet(const Para* q, SlaterDet* Q);
void AlphatoSlaterDet(const Para* q, SlaterDet* Q);
void AlphaprojectgradSlaterDet(const Para* q, 
			       const gradSlaterDet* dQ, double* dq); 


#endif
