/**

  \file ParameterizationClusterFMD.h

  Parametrization of SlaterDet.

  FMD Clusters and FMD nucleons


  (c) 2004 Thomas Neff

*/


#ifndef _PARAMETERIZATIONCLUSTERFMD_H
#define _PARAMETERIZATIONCLUSTERFMD_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationClusterFMD;


int ClusterFMDread(FILE* fp, Para* q);
int ClusterFMDwrite(FILE* fp, const Para* q);

void ClusterFMDclone(const Para* q, Para* qp);

void ClusterFMDinitSlaterDet(const Para* q, SlaterDet* Q);
void ClusterFMDtoSlaterDet(const Para* q, SlaterDet* Q);
void ClusterFMDprojectgradSlaterDet(const Para* q, 
				    const gradSlaterDet* dQ, double* dq); 

#endif
