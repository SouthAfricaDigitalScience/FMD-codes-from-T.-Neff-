/**

  \file MinimizerDONLP2proj.h

  using the minimization routine donlp2


  (c) 2003-2006 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2PROJ_H
#define _MINIMIZERDONLP2PROJ_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/Projection.h"


void initWork(const SlaterDet* Q);

#ifdef MPI
void calcprojectedHamiltonianmpi(const SlaterDet* Q, 
				 const Interaction* Int,
				 int j, int k, int par, 
				 angintegrationpara* projpar,
				 double* eintr, double* eproj);
#endif

void calcprojectedHamiltonian(const SlaterDet* Q, 
			      const Interaction* Int,
			      int j, int k, int par, 
			      angintegrationpara* projpar,
			      double* eintr, double* eproj);


void MinimizeDONLP2proj(const Interaction* Int, int j, int k, int par,
			angintegrationpara* projpar, 
			const Constraint* Const, int nconst,
			Parameterization* P, Para* q,
			int maxsteps, int log, const char* logfile);

#endif
