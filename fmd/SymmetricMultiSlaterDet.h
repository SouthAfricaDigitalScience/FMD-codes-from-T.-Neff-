/**

   \file SymmetricMultiSlaterDet.h

   implement Multiconfig state with Symmetry properties

   (c) 2006 Thomas Neff

*/


#ifndef _SYMMETRICMULTISLATERDET_H
#define _SYMMETRICMULTISLATERDET_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"


extern MultiSlaterDet SymmetricMultiSlaterDet;

Symmetry SymmetricMultiSlaterDetSymmetry(const MultiSlaterDet* MB, int iM);
complex double SymmetricMultiSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i);
void SymmetricMultiSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q);
int SymmetricMultiSlaterDetRead(FILE* fp, MultiSlaterDet* MB);


#endif
