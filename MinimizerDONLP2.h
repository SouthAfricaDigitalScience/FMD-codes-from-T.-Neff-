/**

  \file MinimizerDONLP2.h

  using the minimization routine donlp2


  (c) 2003 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2_H
#define _MINIMIZERDONLP2_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "numerics/donlp2.h"


void MinimizeDONLP2(const Interaction* Int,
		    const Constraint* Const, int nconst,
		    Parameterization* P, Para* q,
		    int maxsteps, int log, const char* logfile);

#endif
