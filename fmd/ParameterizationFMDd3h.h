/**

  \file ParameterizationFMDd3h/h

  Parametrization of SlaterDet with d3h symmetry


  (c) 2009 Thomas Neff

*/


#ifndef _PARAMETERIZATIONFMDD3H_H
#define _PARAMETERIZATIONFMDD3H_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationFMDd3h;


int FMDd3hread(FILE* fp, Para* q);
int FMDd3hwrite(FILE* fp, const Para* q);

void FMDd3hclone(const Para* q, Para* qp);

void FMDd3hinitSlaterDet(const Para* q, SlaterDet* Q);
void FMDd3htoSlaterDet(const Para* q, SlaterDet* Q);
void FMDd3hprojectgradSlaterDet(const Para* q, 
                                const gradSlaterDet* dQ, double* dq); 

#endif
