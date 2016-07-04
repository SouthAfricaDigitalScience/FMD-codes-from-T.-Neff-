/**

  \file ProjectionMultimpi.h

  angular momentum projection of matrix elements
  MPI version


  (c) 2003, 2004, 2005 Thomas Neff

*/


#ifndef _PROJECTIONMULTIMPI_H
#define _PROJECTIONMULTIMPI_H


#include "fmd/SlaterDet.h"
#include "fmd/MultiSlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"


void calcprojectedMultiMBMEmpi(const Projection* P, const ManyBodyOperator* Op,
			       const MultiSlaterDet* MBA, 
			       const MultiSlaterDet* MBB,
			       void** mbme);

#endif
