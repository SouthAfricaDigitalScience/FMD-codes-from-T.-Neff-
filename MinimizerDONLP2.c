/**

  \file MinimizerDONLP2.c

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

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/Constraint.h"
#include "numerics/donlp2.h"

#ifdef MPI
#include "fmdmpi/Hamiltonianmpi.h"
#include "fmdmpi/gradHamiltonianmpi.h"
#else
#include "fmd/Hamiltonian.h"
#include "fmd/gradHamiltonian.h"
#endif

#include "MinimizerDONLP2.h"

#include "numerics/fortranc.h"
#include "misc/physics.h"
#include "misc/utils.h"


// signal handling with MPI jobs is too complicated

#ifndef MPI
jmp_buf env;

volatile sig_atomic_t sigterminate;

void catchterminate(int sig)
{
  sigterminate = 1;
  fprintf(stderr, "... catched signal - terminate minimization\n");
  signal(sig, catchterminate);
}
#endif


// parameters and working space 

static struct {
  int log;
  const char* logfile;
  int maxsteps;
  const Interaction* Int;
  const Constraint* Const;
  int nconst;
  Para* q;
  Parameterization* P;
  SlaterDet* Q;
  SlaterDetAux* X;
  gradSlaterDetAux* dX;
  gradSlaterDet* dH;
} Min;



static void copyxtopara(const double* x, Para* q)
{
  int i;
  for (i=0; i<q->n; i++)
    q->x[i] = x[i];
}


void MinimizeDONLP2(const Interaction* Int, 
		    const Constraint* Const, int nconst,
		    Parameterization* P, Para* q,
		    int maxsteps, int log, const char* logfile)
{
#ifndef MPI
  // handler for INT and TERM signals
  signal(SIGINT, catchterminate);
  signal(SIGTERM, catchterminate);
#endif

  Min.log = log;
  Min.logfile = logfile;
  Min.maxsteps = maxsteps;
  Min.Int = Int;
  Min.Const = Const;
  Min.nconst = nconst;
  Min.q = q;
  Min.P = P;
  Min.Q = (SlaterDet*) malloc(sizeof(SlaterDet));
  Min.X = (SlaterDetAux*) malloc(sizeof(SlaterDetAux));
  Min.dX = (gradSlaterDetAux*) malloc(sizeof(gradSlaterDetAux));
  Min.dH = (gradSlaterDet*) malloc(sizeof(gradSlaterDet));

  Min.P->ParainitSlaterDet(q, Min.Q);
  initSlaterDetAux(Min.Q, Min.X);
  initgradSlaterDetAux(Min.Q, Min.dX);
  initgradSlaterDet(Min.Q, Min.dH);

#ifndef MPI
  // regular program execution or catched signal ?
  if (!setjmp(env))
#endif
    FORTRAN(donlp2)();

  copyxtopara(FORTRAN(o8xdat).x, Min.q);
}


void FORTRAN(ef)(const double* x, double* fx)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAux(Min.Q, Min.X);
#ifdef MPI
  calcHamiltonianmpi(Min.Int, Min.Q, Min.X, fx);
#else
  calcHamiltonian(Min.Int, Min.Q, Min.X, fx);
#endif
  FORTRAN(o8cnt).icf++;
  // fprintf(stderr, "ef:     \tE = %8.3f MeV\n", hbc*(*fx));
}


void FORTRAN(egradf)(const double* x, double* gradf)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAux(Min.Q, Min.X);
  calcgradSlaterDetAux(Min.Q, Min.X, Min.dX);
#ifdef MPI
  calcgradHamiltonianmpi(Min.Int, Min.Q, Min.X, Min.dX, Min.dH);
#else
  calcgradHamiltonian(Min.Int, Min.Q, Min.X, Min.dX, Min.dH);
#endif
  Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradf);
  FORTRAN(o8cnt).icgf++;
  fprintf(stderr, "step %3d: \tE = %8.3f MeV\n", 
	  FORTRAN(o8cnt).icgf, hbc*creal(Min.dH->val));

  if (Min.log && !(FORTRAN(o8cnt).icgf % Min.log)) {
    char logfname[255];
    sprintf(logfname, "%s.min.%03d.log", Min.logfile, FORTRAN(o8cnt).icgf);
    FILE* logfp;
    if (!(logfp = fopen(logfname, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", logfname);
    } else {
      fprintinfo(logfp);
      fprintf(logfp, "# step %3d: \tE = %8.3f MeV\n", 
	      FORTRAN(o8cnt).icgf, hbc*creal(Min.dH->val));
      fprintf(logfp, "\n# Parameterization\n");
      fprintf(logfp, "<Parameterization %s>\n", Min.P->name);
      Min.P->Parawrite(logfp, Min.q);
      fprintf(logfp, "\n# SlaterDet\n");
      writeSlaterDet(logfp, Min.Q);
      fclose(logfp);
    }
  }   

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


void FORTRAN(eh)(int* i, const double* x, double* hxi)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAux(Min.Q, Min.X);
  
  Min.Const[*i-1].me(Min.Q, Min.X, hxi);

  fprintf(stderr, "\t\t\tme   %4s = %8.3f\n", 
  	  Min.Const[*i-1].label, Min.Const[*i-1].output(*hxi));

  *hxi -= Min.Const[*i-1].val;
}


void FORTRAN(egradh)(int* i, const double* x, double* gradhi)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif
  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);
  calcSlaterDetAux(Min.Q, Min.X);
  calcgradSlaterDetAux(Min.Q, Min.X, Min.dX);
  zerogradSlaterDet(Min.dH);
  Min.Const[*i-1].gradme(Min.Q, Min.X, Min.dX, Min.dH);
  Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradhi);
  fprintf(stderr, "\t\t\tgrad %4s = %8.3f\n", 
	  Min.Const[*i-1].label, Min.Const[*i-1].output(Min.dH->val));
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
  FORTRAN(o8par).tau0 = 10.0;  // 20.0
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
