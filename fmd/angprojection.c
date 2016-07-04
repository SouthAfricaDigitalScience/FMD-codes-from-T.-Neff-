/** 

   \file angprojection.c

   Gridpoints used for angular momentum projection

   (c) 2004, 2005 Thomas Neff

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "SlaterDet.h"
#include "Projection.h"
#include "Symmetry.h"

#include "numerics/zcw.h"
#include "numerics/gaussquad.h"


#define SQR(x) ((x)*(x))


int initAngintegration(angintegrationpara* par, const char* projpar)
{
  char projparcpy[strlen(projpar)];
  strcpy(projparcpy, projpar);
  char* c = projparcpy;
	
  c = strtok(c, "-");

  if (!strncmp(c, "ang", 3)) {
    c=strtok(NULL, "-");
    if (!strncmp(c, "none", 4)) {
      par->type = AngNone;
      par->n = 1;
    } else if (!strncmp(c, "zcw", 3)) {
      par->type = AngZCW;
      c=strtok(NULL, "-");
      par->zcw.idx = atoi(c);
      par->n = nangles3(par->zcw.idx);
    } else {
      par->type = AngProd;
      int nbeta = atoi(c);
      c=strtok(NULL, "-");
      int nazimuth = atoi(c);

      par->n = nbeta*SQR(nazimuth);

      par->prod.nbeta = nbeta;
      par->prod.beta = malloc(nbeta*sizeof(double));
      par->prod.wbeta = malloc(nbeta*sizeof(double));
      GaussLegendrePoints(nbeta, -1.0, 1.0, par->prod.beta, par->prod.wbeta);
    
      par->prod.nalpha = nazimuth;
      par->prod.alpha = malloc(nazimuth*sizeof(double));
      par->prod.walpha = malloc(nazimuth*sizeof(double));
      double shift = 0.25;
      ShiftedPeriodicTrapezoidalPoints(nazimuth, 0, 2*M_PI, shift,
				       par->prod.alpha, par->prod.walpha);

      par->prod.ngamma = nazimuth;
      par->prod.gamma = par->prod.alpha;
      par->prod.wgamma = par->prod.walpha;

    }
  }

  return 0;
}


void freeAngintegration(angintegrationpara* par)
{
  if(par->type == AngProd) {
    free(par->prod.beta);
    free(par->prod.wbeta);
    free(par->prod.alpha);
    free(par->prod.walpha);
  }
}


// workspace
static SlaterDet Qpp;
static SlaterDetAux Xpp;

double _estimateangkappa(const SlaterDet* Q)
{
  const double cosb0 = 0.5;

  double ovl0;
  complex double ovlcosb0;
   
  // do we have to (re)initialize workspace first ?
  if (Qpp.A < Q->A) {
    initSlaterDet(Q, &Qpp);
    initSlaterDetAux(Q, &Xpp);
  }

  copySlaterDet(Q, &Qpp);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);
  ovl0 = Xpp.ovlap;

  rotateSlaterDet(&Qpp, 0.0,acos(cosb0),0.0);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);
  ovlcosb0 = Xpp.ovlap;

  return (log(cabs(ovlcosb0/ovl0))/log(cosb0));
}
  

// kappacrit
static double kappacrit = 25.0;

void _setangkappacrit(double kappa)
{
  kappacrit = kappa;
}

double _getangkappacrit()
{
  return (kappacrit);
}


void _initangintegration(const Projection* P, 
                         double kappa,
                         Symmetry S, Symmetry Sp,
                         angintegrationpara* par)
{
  par->type = P->ang;

  if (P->ang == AngNone) {
    par->n = 1;
  }

  else if (P->ang == AngZCW) {
    par->zcw.idx = P->angzcw.idx;
    par->n = nangles3(P->angzcw.idx);
  }

  else if (P->ang == AngProd || P->ang == AngProdA) {
    int nbeta = P->angprod.nbeta;
    int nazim = P->angprod.nazimuth;

    // use Gauss-Legendre integration for cos(beta)
    // if overlap is very strongly peaked, use Gauss-Exponential 

    if (hasSymmetry(S, spherical))
      par->prod.nbeta = 1;
    else	
      par->prod.nbeta = nbeta;    

    par->prod.beta = malloc(par->prod.nbeta*sizeof(double));
    par->prod.wbeta = malloc(par->prod.nbeta*sizeof(double));

    if (kappa < kappacrit)
      GaussLegendrePoints(par->prod.nbeta, -1.0, 1.0, par->prod.beta, par->prod.wbeta);
    else
      GaussExponentialPoints(par->prod.nbeta, -1.0, 1.0, kappa, par->prod.beta, par->prod.wbeta);

    // azimuthal integration can be skipped if symmetry allows
    // trapezoidal rule with shifted grid points
    // special handling for rotatez2, rotatez3 ?
    double shift = 0.25;

    if (hasAxialSymmetry(S))
      par->prod.nalpha = 1;
    else
      par->prod.nalpha = nazim;

    par->prod.alpha = malloc(par->prod.nalpha*sizeof(double));
    par->prod.walpha = malloc(par->prod.nalpha*sizeof(double));

    ShiftedPeriodicTrapezoidalPoints(par->prod.nalpha, 0, 2*M_PI, shift,
				     par->prod.alpha, par->prod.walpha);

    if (hasAxialSymmetry(Sp))
      par->prod.ngamma = 1;
    else
      par->prod.ngamma = nazim;

    par->prod.gamma = malloc(par->prod.ngamma*sizeof(double));
    par->prod.wgamma = malloc(par->prod.ngamma*sizeof(double));

    ShiftedPeriodicTrapezoidalPoints(par->prod.ngamma, 0, 2*M_PI, shift,
				     par->prod.gamma, par->prod.wgamma);

    par->n = par->prod.nbeta* par->prod.nalpha* par->prod.ngamma;
  }
}


static inline double dmin(double a, double b)
{
  return (a < b ? a : b);
}


// What's the proper criterium for switchung to Gauss-Exponential ?
// has been overlap <Q|R(beta)|Qp>
// this does not work for Gamow-Teller and Spectroscopic Amplitudes
// now use min(<Q|R(beta)|Q>, <Qp|R(beta)|Qp>)

void initangintegration(const Projection* P, 
			const SlaterDet* Q, const SlaterDet* Qp,
			Symmetry S, Symmetry Sp,
			angintegrationpara* par)
{
  double kappa = 0.0; 
  if (P->ang == AngProdA) {
    kappa = dmin(_estimateangkappa(Q), _estimateangkappa(Qp));
    fprintf(stderr, "kappa: %8.3f\n", kappa);
  }

  _initangintegration(P, kappa,  S, Sp, par);
}


void reinitangintegration(const SlaterDet* Q, const SlaterDet* Qp,
                          Symmetry S, Symmetry Sp,
                          angintegrationpara* par)
{
  // we only have to modify something for adaptive integration

  // use Gauss-Legendre integration for cos(beta)
  // if overlap is very strongly peaked, use Gauss-Exponential 

  double kappa = 0.0;

  if (par->type == AngProdA) {

    kappa = dmin(_estimateangkappa(Q), _estimateangkappa(Qp));
    fprintf(stderr, "estimated kappa: %8.3f\n", kappa);

    if (kappa < kappacrit)
      GaussLegendrePoints(par->prod.nbeta, -1.0, 1.0, par->prod.beta, par->prod.wbeta);
    else
      GaussExponentialPoints(par->prod.nbeta, -1.0, 1.0, kappa, par->prod.beta, par->prod.wbeta);

  }
}


void getangintegrationpoint(int i, const angintegrationpara* par,
			    double* alpha, double* beta, double* gamma, 
			    double* w)
{
  assert(i < par->n);

  if (par->type == AngNone) {
    *alpha = 0.0; *beta = 0.0; *gamma = 0.0;
    *w = 8*SQR(M_PI);
  }

  else if (par->type == AngZCW) {
    getangles3(par->zcw.idx, i, alpha, beta, gamma);
    *w = 8*SQR(M_PI)/par->n;
  }

  else if (par->type == AngProd || par->type == AngProdA) {
    int ibeta, ialpha, igamma;
    
    ibeta = i % par->prod.nbeta;
    ialpha = (i/par->prod.nbeta) % par->prod.nalpha;
    igamma = i/(par->prod.nbeta*par->prod.nalpha);

    *beta = acos(par->prod.beta[ibeta]);
    *alpha = par->prod.alpha[ialpha];
    *gamma = par->prod.gamma[igamma];
    *w = par->prod.wbeta[ibeta]*
      par->prod.walpha[ialpha]*par->prod.wgamma[igamma];
  }
    
  // fprintf(stderr, "ang integration point: (%8.5f, %8.5f, %8.5f) weight %8.5f\n", 
  //	     *alpha, *beta, *gamma, *w);
}
