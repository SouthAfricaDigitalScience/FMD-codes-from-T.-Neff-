/**

  \file MinimizerDONLP2p.h

  using the minimization routine donlp2


  (c) 2003 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2P_H
#define _MINIMIZERDONLP2P_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"


void MinimizeDONLP2p(const Interaction* Int, int par, 
		     const Constraint* Const, int nconst,
		     Parameterization* P, Para* q,
		     int maxsteps, int log, const char* logfile);

#endif
