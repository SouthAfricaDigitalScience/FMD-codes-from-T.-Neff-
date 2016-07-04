/**

  \file MinimizerDONLP2par.c

  using the minimization routine donlp2

  this is ugly as hell, but it works !


  (c) 2003 Thomas Neff

*/


// we need some information about the FMD parameterization
// regarding the number of parameters per Gaussian and about
// the position of the parameter a

#define NGAUSS 12
#define APOS 4

#include <signal.h>
#include <stdlib.h>
#include <setjmp.h>
#include <string.h>
#include <math.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "fmd/ConstraintCM.h"
#include "numerics/donlp2.h"

#ifdef MPI
#include "fmdmpi/Hamiltonianmpi.h"
#include "fmdmpi/gradHamiltonianmpi.h"
#else
#include "fmd/Hamiltonian.h"
#include "fmd/gradHamiltonian.h"
#endif

#include "MinimizerDONLP2par.h"

#include "numerics/fortranc.h"
#include "misc/physics.h"


jmp_buf env;

volatile sig_atomic_t sigterminate;

void catchterminate(int sig)
{
  sigterminate = 1;
  fprintf(stderr, "... catched signal - terminate minimization\n");
  signal(sig, catchterminate);
}


// parameters and working space 

static struct {
  int maxsteps;
  const Interaction* Int;
  int par;
  const Constraintod* Const;
  int nconst;
  Para* q;
  Parameterization* P;
  SlaterDet* Q;
  SlaterDet* Qp;
  SlaterDetAux* X;
  SlaterDetAux* Xp;
  gradSlaterDetAux* dX;
  gradSlaterDet* dH, *dHd, *dHp, *dNd, *dNp;
} Min;



static void copyxtopara(const double* x, Para* q)
{
  int i;
  for (i=0; i<q->n; i++)
    q->x[i] = x[i];
}


void MinimizeDONLP2par(const Interaction* Int, int par,
		       const Constraintod* Const, int nconst,
		       Parameterization* P, Para* q,
		       int maxsteps)
{
  // handler for INT and TERM signals
  signal(SIGINT, catchterminate);
  signal(SIGTERM, catchterminate);

  Min.maxsteps = maxsteps;
  Min.Int = Int;
  Min.par = par;
  Min.Const = Const;
  Min.nconst = nconst;
  Min.q = q;
  Min.P = P;
  Min.Q = (SlaterDet*) malloc(sizeof(SlaterDet));
  Min.Qp = (SlaterDet*) malloc(sizeof(SlaterDet));
  Min.X = (SlaterDetAux*) malloc(sizeof(SlaterDetAux));
  Min.Xp = (SlaterDetAux*) malloc(sizeof(SlaterDetAux));
  Min.dX = (gradSlaterDetAux*) malloc(sizeof(gradSlaterDetAux));
  Min.dH = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));
  Min.dHd = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));
  Min.dHp = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));
  Min.dNd = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));
  Min.dNp = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));

  Min.P->ParainitSlaterDet(q, Min.Q);
  Min.P->ParainitSlaterDet(q, Min.Qp);
  initSlaterDetAux(Min.Q, Min.X);
  initSlaterDetAux(Min.Q, Min.Xp);
  initgradSlaterDetAux(Min.Q, Min.dX);
  initgradSlaterDet(Min.Q, Min.dH);
  initgradSlaterDet(Min.Q, Min.dHd);
  initgradSlaterDet(Min.Q, Min.dHp);
  initgradSlaterDet(Min.Q, Min.dNd);
  initgradSlaterDet(Min.Q, Min.dNp);

  // regular program execution or catched signal ?
  if (!setjmp(env))
    FORTRAN(donlp2)();

  copyxtopara(FORTRAN(o8xdat).x, Min.q);
}


void FORTRAN(ef)(const double* x, double* fx)
{
  complex double hd, hp, nd, np;

  if (sigterminate)
    longjmp(env, 1);
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAuxod(Min.Q, Min.Q, Min.X);
  nd = Min.X->ovlap;
#ifdef MPI
  calcHamiltonianodmpi(Min.Int, Min.Q, Min.Q, Min.X, &hd);
#else
  calcHamiltonianod(Min.Int, Min.Q, Min.Q, Min.X, &hd);
#endif
  copySlaterDet(Min.Q, Min.Qp);
  invertSlaterDet(Min.Qp);
  calcSlaterDetAuxod(Min.Q, Min.Qp, Min.X);
  np = Min.X->ovlap;
#ifdef MPI
  calcHamiltonianodmpi(Min.Int, Min.Q, Min.Qp, Min.X, &hp);
#else
  calcHamiltonianod(Min.Int, Min.Q, Min.Qp, Min.X, &hp);
#endif
  *fx = (hd+Min.par*hp)/(nd+Min.par*np);

  FORTRAN(o8cnt).icf++;
}


void FORTRAN(egradf)(const double* x, double* gradf)
{
  complex double hd, hp, nd, np;
  complex double H, N;

  if (sigterminate)
    longjmp(env, 1);
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAuxod(Min.Q, Min.Q, Min.X);
  calcgradSlaterDetAuxod(Min.Q, Min.Q, Min.X, Min.dX);
  calcgradOvlapod(Min.Q, Min.Q, Min.X, Min.dX, Min.dNd);
#ifdef MPI
  calcgradHamiltonianodmpi(Min.Int, Min.Q, Min.Q, Min.X, Min.dX, Min.dHd);
#else
  calcgradHamiltonianod(Min.Int, Min.Q, Min.Q, Min.X, Min.dX, Min.dHd);
#endif
  hd = Min.dHd->val;
  nd = Min.X->ovlap;

  copySlaterDet(Min.Q, Min.Qp);
  invertSlaterDet(Min.Qp);
  calcSlaterDetAuxod(Min.Q, Min.Qp, Min.X);
  calcgradSlaterDetAuxod(Min.Q, Min.Qp, Min.X, Min.dX);
  calcgradOvlapod(Min.Q, Min.Qp, Min.X, Min.dX, Min.dNp);
#ifdef MPI
  calcgradHamiltonianodmpi(Min.Int, Min.Q, Min.Qp, Min.X, Min.dX, Min.dHp);
#else
  calcgradHamiltonianod(Min.Int, Min.Q, Min.Qp, Min.X, Min.dX, Min.dHp);
#endif
  hp = Min.dHp->val;
  np = Min.X->ovlap;

  H = hd+Min.par*hp;
  N = nd+Min.par*np;
    
  zerogradSlaterDet(Min.dH);
  addmulttogradSlaterDet(Min.dH, Min.dHd, 1.0/N);
  addmulttogradSlaterDet(Min.dH, Min.dHp, Min.par* 1.0/N);
  addmulttogradSlaterDet(Min.dH, Min.dNd, -H/(N*N));
  addmulttogradSlaterDet(Min.dH, Min.dNp, -Min.par*H/(N*N));
  
  Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradf);
  FORTRAN(o8cnt).icgf++;
  fprintf(stderr, "step %3d: \tE = %8.3f (%8.3f) MeV\n", 
	  FORTRAN(o8cnt).icgf, hbc*creal(H/N), hbc*creal(hd/nd));
}


void FORTRAN(eg)(int* i, const double* x, double* gxi)
{
  if (!strcmp(Min.P->name, "FMD")) {
    *gxi = x[(*i-1)*NGAUSS+APOS];
  }
}


void FORTRAN(egradg)(int* i, const double* x, double* gradgi)
{
}


#define SQR(x) ((x)*(x))

void FORTRAN(eh)(int* i, const double* x, double* hxi)
{
  if (sigterminate)
    longjmp(env, 1);

  // First Constraint is CM constraint
  if (*i == 1) {

    Constraint* cnst = &ConstraintCM;

    copyxtopara(x, Min.q);
    Min.P->ParatoSlaterDet(Min.q, Min.Q);
    calcSlaterDetAux(Min.Q, Min.X);
  
    cnst->me(Min.Q, Min.X, hxi);

  } else {    

    Constraintod* cnst = &Min.Const[*i-1];

    complex double nd, np;
    complex double med[cnst->dim], mep[cnst->dim];
    double expme[cnst->dim];
    int j;

    copyxtopara(x, Min.q);
    Min.P->ParatoSlaterDet(Min.q, Min.Q);
    calcSlaterDetAuxod(Min.Q, Min.Q, Min.X);
    nd = Min.X->ovlap;
    copySlaterDet(Min.Q, Min.Qp);
    invertSlaterDet(Min.Qp);
    calcSlaterDetAuxod(Min.Q, Min.Qp, Min.Xp);
    np = Min.Xp->ovlap;
  
    cnst->me(Min.Q, Min.Q, Min.X, med);
    cnst->me(Min.Q, Min.Qp, Min.Xp, mep);

    for (j=0; j<cnst->dim; j++)
      expme[j] = (med[j] + Min.par* mep[j])/(nd + Min.par* np);

    if (cnst->dim==1)
      *hxi = expme[0];
    else {
      double hctcr = 0.0;
      for (j=0; j<cnst->dim; j++)
	hctcr += SQR(expme[j]);
      *hxi = sqrt(hctcr);
    }

    *hxi -= cnst->val;
  }
}


void FORTRAN(egradh)(int* i, const double* x, double* gradhi)
{
  if (sigterminate)
    longjmp(env, 1);

  // First constraint is CM constraint
  if (*i == 1) {
    Constraint* cnst = &ConstraintCM;

    copyxtopara(x, Min.q);

    Min.P->ParatoSlaterDet(Min.q, Min.Q);
    calcSlaterDetAux(Min.Q, Min.X);
    calcgradSlaterDetAux(Min.Q, Min.X, Min.dX);
    zerogradSlaterDet(Min.dH);
    cnst->gradme(Min.Q, Min.X, Min.dX, Min.dH);
    Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradhi);
    fprintf(stderr, "\t\t\t\t\t%4s = %8.3f\n", 
	    cnst->label, cnst->output(Min.dH->val));
  } else {
    Constraintod* cnst = &Min.Const[*i-1];

    complex double nd, np, hd, hp;
    complex double N, H;
    double expme[cnst->dim];

    copyxtopara(x, Min.q);
    Min.P->ParatoSlaterDet(Min.q, Min.Q);
    copySlaterDet(Min.Q, Min.Qp);
    invertSlaterDet(Min.Qp);

    calcSlaterDetAuxod(Min.Q, Min.Q, Min.X);
    nd = Min.X->ovlap;
    calcSlaterDetAuxod(Min.Q, Min.Qp, Min.Xp);
    np = Min.Xp->ovlap;

    N = nd+Min.par*np;

    if (cnst->gradneedsme) {
      complex double med[cnst->dim], mep[cnst->dim];

      cnst->me(Min.Q, Min.Q, Min.X, med);
      cnst->me(Min.Q, Min.Qp, Min.Xp, mep);

      int j;
      for (j=0; j<cnst->dim; j++)
	expme[j] = (med[j] + Min.par* mep[j])/N;
    }

    calcgradSlaterDetAuxod(Min.Q, Min.Q, Min.X, Min.dX);
    calcgradOvlapod(Min.Q, Min.Q, Min.X, Min.dX, Min.dNd);
    zerogradSlaterDet(Min.dHd);
    cnst->gradme(Min.Q, Min.Q, Min.X, Min.dX, expme, Min.dHd);

    hd = Min.dHd->val;

    calcgradSlaterDetAuxod(Min.Q, Min.Qp, Min.Xp, Min.dX);
    calcgradOvlapod(Min.Q, Min.Qp, Min.Xp, Min.dX, Min.dNp);
    zerogradSlaterDet(Min.dHp);
    cnst->gradme(Min.Q, Min.Qp, Min.Xp, Min.dX, expme, Min.dHp);

    hp = Min.dHp->val;

    H = hd+Min.par*hp;
    
    zerogradSlaterDet(Min.dH);
    addmulttogradSlaterDet(Min.dH, Min.dHd, 1.0/N);
    addmulttogradSlaterDet(Min.dH, Min.dHp, Min.par* 1.0/N);
    addmulttogradSlaterDet(Min.dH, Min.dNd, -H/(N*N));
    addmulttogradSlaterDet(Min.dH, Min.dNp, -Min.par*H/(N*N));
  
    Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradhi);
    fprintf(stderr, "\t\t\t\t\t%4s = %8.3f (%8.3f)\n",
	    cnst->label, cnst->output(creal(H/N)), cnst->output(creal(hd/nd)));

  }
}


void FORTRAN(setup0)(void)
{
  int i;

  FORTRAN(o8dim).n = Min.q->n;
  FORTRAN(o8dim).nh = Min.nconst;
  // in case of FMD parametrization use bound constraints for Re{a}
  if (!strcmp(Min.P->name, "FMD")) {
    FORTRAN(o8dim).ng = Min.q->ngauss;
    for (i=0; i<Min.nconst+1; i++) {
      FORTRAN(o8gri).gunit[i][0] = -1;
      FORTRAN(o8gri).gunit[i][1] = 0;
      FORTRAN(o8gri).gunit[i][2] = 0;
    }
    for (i=0; i<Min.q->ngauss; i++) {
      FORTRAN(o8gri).gunit[i+Min.nconst+1][0] = 1;
      FORTRAN(o8gri).gunit[i+Min.nconst+1][1] = i*NGAUSS+APOS+1;
      FORTRAN(o8gri).gunit[i+Min.nconst+1][2] = 1;
    }
    for (i=0; i<Min.q->ngauss; i++) {
      if (Min.q->x[i*NGAUSS+APOS] < 0)
	Min.q->x[i*NGAUSS+APOS] = 1.0+1.0*rand()/RAND_MAX;
    }
  } else
    FORTRAN(o8dim).ng = 0;
  FORTRAN(o8stpa).analyt = 1;
  FORTRAN(o8stpa).cold = 1;
  FORTRAN(o8par).tau0 = 20.0;
  FORTRAN(o8par).del0 = 0.0;
  FORTRAN(o8stpa).silent = 1;

  for (i=0; i<Min.q->n; i++)
    FORTRAN(o8xdat).x[i] = Min.q->x[i];
}


void FORTRAN(setup)(void)
{
  FORTRAN(o8stpa).te0 = 0;
  FORTRAN(o8stpa).te1 = 0;
  FORTRAN(o8par).iterma = Min.maxsteps-1;
}
