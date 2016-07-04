/**

  \file Projectionmpi.h

  angular momentum projection of matrix elements
  MPI version


  (c) 2003, 2004, 2005 Thomas Neff

*/


#ifndef _PROJECTIONMPI_H
#define _PROJECTIONMPI_H


#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"


void calcprojectedMBMEmpi(const Projection* P, const ManyBodyOperator* Op,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  Symmetry S, Symmetry Sp,
			  void* mbme);

void calcprojectedMBMEsmpi(const Projection* P, const ManyBodyOperators* Ops,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   Symmetry S, Symmetry Sp,
			   void* mbme);


#endif
