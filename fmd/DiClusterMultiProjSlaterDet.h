/**

   \file DiClusterMultiProjSlaterDet.h

   set of MultiSlaterDet built out of two Clusters

   each cluster is a multiconfig state projected on a set of (J, Jz) values

   
   (c) 2006 Thomas Neff

*/


#ifndef _DICLUSTERMULTIPROJSLATERDET_H
#define _DICLUSTERMULTIPROJSLATERDET_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet DiClusterMultiProjSlaterDet;

Symmetry DiClusterMultiProjSlaterDetSymmetry(const MultiSlaterDet* MSD, int iM);
complex double DiClusterMultiProjSlaterDetWeight(const MultiSlaterDet* MSD, int iM, int i);
void DiClusterMultiProjSlaterDetGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q);
int DiClusterMultiProjSlaterDetRead(FILE* fp, MultiSlaterDet* MSD);


#endif
