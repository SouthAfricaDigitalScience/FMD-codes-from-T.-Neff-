/**

   \file DiClusterMulticonfig.h

   set of MultiSlaterDet built out of two Clusters

   each cluster is a multiconfig state projected on a set of (J, Jz) values

   
   (c) 2009 Thomas Neff

*/


#ifndef _DICLUSTERMULTICONFIG_H
#define _DICLUSTERMULTICONFIG_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet DiClusterMulticonfig;

Symmetry DiClusterMulticonfigSymmetry(const MultiSlaterDet* MSD, int iM);
complex double DiClusterMulticonfigWeight(const MultiSlaterDet* MSD, int iM, int i);
void DiClusterMulticonfigGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q);
int DiClusterMulticonfigRead(FILE* fp, MultiSlaterDet* MSD);


#endif
