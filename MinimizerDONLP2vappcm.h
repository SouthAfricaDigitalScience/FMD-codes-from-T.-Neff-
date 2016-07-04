/**

  \file MinimizerDONLP2vappcm.h

  using the minimization routine donlp2


  (c) 2003-2011 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2VAPPCM_H
#define _MINIMIZERDONLP2VAPPCM_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/Projection.h"


void initWork(int j, int par,
              angintegrationpara* angpara,
	      cmintegrationpara* cmpara,
              const Interaction* Int,
              const SlaterDet* Q);


void calcprojectedHamiltonian(const SlaterDet* Q, 
                              const Interaction* Int,
                              int j, int par, 
                              int ival, 
                              angintegrationpara* angpara,
			      cmintegrationpara* cmpara,
                              double* eintr, double* eproj);


void MinimizeDONLP2vapp(const Interaction* Int, 
			int j, int par, int ival,
			double threshkmix, double minnormkmix,
			double alpha,
			angintegrationpara* angpara,
			cmintegrationpara* cmpara,
			const Constraint* Const, int nconst,
			Parameterization* P, Para* q,
			int maxsteps, int log, const char* logfile);

#endif
