/**

  \file MinimizerDONLP2vapp.h

  using the minimization routine donlp2


  (c) 2003-2006 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2VAPP_H
#define _MINIMIZERDONLP2VAPP_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/Projection.h"


void initWork(int j, int par,
              angintegrationpara* angpara,
              const Interaction* Int,
              const SlaterDet* Q);


void calcprojectedHamiltonian(const SlaterDet* Q, 
                              const Interaction* Int,
                              int j, int par, 
                              int ival, 
                              angintegrationpara* angpara,
                              double* eintr, double* eproj);


void MinimizeDONLP2vapp(const Interaction* Int, 
                        int j, int par, int ival,
                        double threshkmix, double minnormkmix,
                        double alpha,
                        angintegrationpara* projpar, 
                        const Constraint* Const, int nconst,
                        Parameterization* P, Para* q,
                        int maxsteps, int log, const char* logfile);

#endif
