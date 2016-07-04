/**

  \file donlp2

  the minimization routine donlp2 by Peter Spellucci


  (c) 2003 Thomas Neff

*/


#ifndef _DONLP2_H
#define _DONLP2_H


#include "fortranc.h"

// donlp2 compiled with 
// NX=1200,NRESM=100,MAXIT=10000,NSTEP=40
// keep in sync with O8PARA.INC !

#define NX 1200
#define NRESM 100
#define MAXIT 10000
#define NSTEP 40

// donlp2's common blocks

extern struct {
  double x[NX];
  double x0[NX];
  double x1[NX];
  double xmin[NX];
  // further entries
} FORTRAN(o8xdat);


extern struct {
  double xst[NX];
} FORTRAN(o8stv);

extern struct {
  int nreset;
  int numsm;
} FORTRAN(o8rst);

extern struct {
  int n;
  int nh;
  int ng;
  int nr;
  int nres;
} FORTRAN(o8dim);

extern struct {
  int val[NRESM+1];
  int gconst[NRESM+1];
  int gunit[NRESM+1][3];
  int LLOW[NX];
  int LUP[NX];
} FORTRAN(o8gri);

extern struct {
  double del0;
  double del01;
  double del;
  double delmin;
  double tau0;
  double tau;
  double smalld;
  double smallw;
  double rho;
  double rho1;
  double eta;
  double ny;
  double epsx;
  double epshi;
  double c1d;
  double scfmax;
  double tauqp;
  double taufac;
  double taumax;
  double updmy0;
  int iterma;
  int ifill1;
} FORTRAN(o8par);

// logicals are ints
extern struct {
  int intakt;
  int inx;
  int std;
  int te0;
  int te1;
  int te2;
  int te3;
  int singul;
  int ident;
  int eqres;
  int silent;
  int analyt;
  int cold;
} FORTRAN(o8stpa);

extern struct {
  int icf;
  int icgf;
  int cfincr;
  int cres[NRESM];
  int cgres[NRESM];
} FORTRAN(o8cnt);

extern struct {
  float optite;
  int itstep;
  int phase;
  float runtim;
  double accinf[32][MAXIT+1];
} FORTRAN(o8itin);

extern struct {
  int ffuerr;
  int cfuerr[NRESM];
} FORTRAN(o8err);


// the minimizer itself

void FORTRAN(donlp2)(void);

// these routines have to be defined for donlp2

void FORTRAN(ef)(const double* x, double* fx);

void FORTRAN(egradf)(const double* x, double* gradf);

void FORTRAN(eg)(int* i, const double* x, double* gxi);

void FORTRAN(egradg)(int* i, const double* x, double* gradgi);

void FORTRAN(eh)(int* i, const double* x, double* hxi);

void FORTRAN(egradh)(int* i, const double* x, double* gradhi);

void FORTRAN(setup0)(void);

void FORTRAN(setup)(void);


#endif
