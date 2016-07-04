/**

   \file DiClusterProjSlaterDet.h

   implement Di-Cluster SlaterDet as MultiSlaterDet

   
   (c) 2005

*/


#ifndef _DICLUSTERPROJSLATERDET_H
#define _DICLUSTERPROJSLATERDET_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet DiClusterProjSlaterDet;

Symmetry DiClusterProjSlaterDetSymmetry(const MultiSlaterDet* MSD, int iM);
complex double DiClusterProjSlaterDetWeight(const MultiSlaterDet* MSD, int iM, int i);
void DiClusterProjSlaterDetGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q);
int DiClusterProjSlaterDetRead(FILE* fp, MultiSlaterDet* MSD);


#endif
