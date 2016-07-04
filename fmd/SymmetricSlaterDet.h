/**

   \file SymmetricSlaterDet.h

   implement SlaterDet with Symmetry properties

   (c) 2005

*/


#ifndef _SYMMETRICSLATERDET_H
#define _SYMMETRICSLATERDET_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet SymmetricSlaterDet;

Symmetry SymmetricSlaterDetSymmetry(const MultiSlaterDet* MB, int iM);
complex double SymmetricSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i);
void SymmetricSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q);
int SymmetricSlaterDetRead(FILE* fp, MultiSlaterDet* MB);


#endif
