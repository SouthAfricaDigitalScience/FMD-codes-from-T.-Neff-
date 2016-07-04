/**

  \file MinimizerDONLP2orthogonalvapp.h

  using the minimization routine donlp2


  (c) 2003-2006 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2ORTHOGONALVAPP_H
#define _MINIMIZERDONLP2ORTHOGONALVAPP_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/Projection.h"


void initWork(int nfix,
              const SlaterDet* Qfix,
              int j, int par,
              angintegrationpara* angpara,
              const Interaction* Int,
              const SlaterDet* Q);


void calcprojectedHamiltonian(int nfix,
                              const SlaterDet* Qfix,
                              const SlaterDet* Q, 
                              const Interaction* Int,
                              int j, int par, 
                              int ifix, 
                              angintegrationpara* angpara,
                              double* eintr, double* emulti);


void MinimizeDONLP2orthogonalvapp(const Interaction* Int, 
				  int j, int par, int ival, double alpha,
				  angintegrationpara* projpar, 
				  const Constraint* Const, int nconst,
				  int ifix, SlaterDet* Qfix,
				  Parameterization* P, Para* q,
				  int maxsteps, int log, const char* logfile);

#endif
