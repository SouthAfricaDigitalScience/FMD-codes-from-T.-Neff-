/**

  \file Overlap.h

  calculate Overlap for Slater determinants


  (c) 2003, 2005 Thomas Neff

*/


#ifndef _OVERLAP_H
#define _OVERLAP_H

#include <stdio.h>

#include "SlaterDet.h"
#include "Projection.h"


extern ManyBodyOperator OpOvlap;


void calcOvlapod(void* dummy,
		 const SlaterDet* Q, const SlaterDet* Qp,
		 const SlaterDetAux* X,
		 complex double* ovl);


void writeprojectedtransitionOvlap(FILE* fp,
                                   const Projection* P,
                                   const complex double**** ovlme,
                                   const Eigenstates* Efin, 
                                   const Eigenstates* Eini);

#endif
