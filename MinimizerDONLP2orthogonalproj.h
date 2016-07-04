/**

  \file MinimizerDONLP2orthogonalproj.h

  using the minimization routine donlp2


  (c) 2006 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2ORTHOGONALPROJ_H
#define _MINIMIZERDONLP2ORTHOGONALPROJ_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/Projection.h"


void initWork(const SlaterDet* Q);


void calcprojectedHamiltonian(const SlaterDet* Q0,
			      const Interaction* Int,
			      int j, int k, int par, 
			      angintegrationpara* projpar,
			      double* h0, double* n0);

#ifdef MPI
void calcprojectedHamiltonianmpi(const SlaterDet* Q0,
				 const Interaction* Int,
				 int j, int k, int par, 
				 angintegrationpara* projpar,
				 double* h0, double* n0);
#endif

void calcorthprojectedHamiltonian(const SlaterDet* Q0, 
				  double h0, double n0,
				  const SlaterDet* Q,
				  const Interaction* Int,
				  int j, int k, int par, 
				  angintegrationpara* projpar,
				  double* eintr, double* eproj);

#ifdef MPI
void calcorthprojectedHamiltonianmpi(const SlaterDet* Q0, 
				     double h0, double n0,
				     const SlaterDet* Q,
				     const Interaction* Int,
				     int j, int k, int par, 
				     angintegrationpara* projpar,
				     double* eintr, double* eproj);
#endif

void MinimizeDONLP2orthproj(const Interaction* Int, int j, int k, int par, 
			    angintegrationpara* projpar,
			    SlaterDet* Q0, double h0, double n0,
			    const Constraint* Const, int nconst,
			    Parameterization* P, Para* q,
			    int maxsteps, int log, const char* logfile);

#endif
