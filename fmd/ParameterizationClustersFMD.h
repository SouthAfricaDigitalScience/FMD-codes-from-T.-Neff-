/**

  \file ParameterizationClustersFMD.h

  Parametrization of SlaterDet.

  FMD Clusters (position, momentum, orientation) and FMD nucleons


  (c) 2006 Thomas Neff

*/


#ifndef _PARAMETERIZATIONCLUSTERSFMD_H
#define _PARAMETERIZATIONCLUSTERSFMD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationClustersFMD;


int ClustersFMDread(FILE* fp, Para* q);
int ClustersFMDwrite(FILE* fp, const Para* q);

void ClustersFMDclone(const Para* q, Para* qp);

void ClustersFMDinitSlaterDet(const Para* q, SlaterDet* Q);
void ClustersFMDtoSlaterDet(const Para* q, SlaterDet* Q);
void ClustersFMDprojectgradSlaterDet(const Para* q, 
				    const gradSlaterDet* dQ, double* dq); 

#endif
