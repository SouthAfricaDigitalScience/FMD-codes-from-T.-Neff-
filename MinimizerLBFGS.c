/*

  Minimizer.c

  (c) 2003 Thomas Neff

*/

#include <stdlib.h>
#include <float.h>

#include "MinimizerLBFGS.h"

#include "numerics/lbfgs.h"


void initMinimizer(Minimizer* mini, int n, double precision)
{
  mini->n = n;;
  mini->m = 6;
  mini->diag = (double*) malloc(mini->n*sizeof(double));
  mini->w = (double*) malloc((mini->n*(2*mini->m+1)+2*mini->m)*sizeof(double));
  mini->diagco = 0;
  mini->iprint[0] = -1;
  mini->iprint[1] = 0;
  mini->eps = precision;
  mini->xtol = DBL_EPSILON;
  mini->iflag = 0;
}


int MinimizerStep(Minimizer* mini, double* x, const double H, const double* dH)
{
  FORTRAN(lbfgs)(&mini->n, &mini->m, x, &H, dH,
		 &mini->diagco, mini->diag, mini->iprint, 
		 &mini->eps, &mini->xtol, mini->w, &mini->iflag);

  return (mini->iflag == 1);
}
