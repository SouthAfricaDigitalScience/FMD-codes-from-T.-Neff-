/**

  \file MinimzerLBFGS.h


  interface to lbfgs minimizer

  (c) 2003 Thomas Neff

*/

#ifndef _MINIMIZERLBFGS_H
#define _MINIMIZERLBFGS_H


typedef struct {
  int n;
  int m;
  int iprint[2];
  int iflag;
  double* diag;
  double* w;
  double eps;
  double xtol;
  unsigned int diagco;
} Minimizer;


void initMinimizer(Minimizer* mini, int n, double precision);

int MinimizerStep(Minimizer* mini, double *x, 
		  const double H, const double* dH);


#endif
