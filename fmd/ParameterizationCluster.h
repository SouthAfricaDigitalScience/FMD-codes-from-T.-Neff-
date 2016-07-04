/**

  \file ParameterizationCluster.h

  Parametrization of SlaterDet.

  FMD Clusters


  (c) 2004 Thomas Neff

*/


#ifndef _PARAMETERIZATIONCLUSTER_H
#define _PARAMETERIZATIONCLUSTER_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationCluster;


int Clusterread(FILE* fp, Para* q);
int Clusterwrite(FILE* fp, const Para* q);

void Clusterclone(const Para* q, Para* qp);

void ClusterinitSlaterDet(const Para* q, SlaterDet* Q);
void ClustertoSlaterDet(const Para* q, SlaterDet* Q);
void ClusterprojectgradSlaterDet(const Para* q, 
				 const gradSlaterDet* dQ, double* dq); 

#endif
