/**

   \file SimpleSlaterDet.h

   implement SimpleSlaterDet as MultiSlaterDet

   (c) 2005

*/


#ifndef _SIMPLESLATERDET_H
#define _SIMPLESLATERDET_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet SimpleSlaterDet;

Symmetry SimpleSlaterDetSymmetry(const MultiSlaterDet* MB, int iM);
complex double SimpleSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i);
void SimpleSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q);
int SimpleSlaterDetRead(FILE* fp, MultiSlaterDet* MB);


#endif
