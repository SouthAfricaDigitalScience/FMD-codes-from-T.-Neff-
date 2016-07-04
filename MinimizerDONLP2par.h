/**

  \file MinimizerDONLP2par.h

  using the minimization routine donlp2


  (c) 2003 Thomas Neff

*/


#ifndef _MINIMIZERDONLP2PAR_H
#define _MINIMIZERDONLP2PAR_H

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"


void MinimizeDONLP2par(const Interaction* Int, int par, 
		       const Constraintod* Const, int nconst,
		       Parameterization* P, Para* q,
		       int maxsteps);

#endif
