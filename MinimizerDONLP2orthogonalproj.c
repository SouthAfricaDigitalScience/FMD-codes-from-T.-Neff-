/**

  \file MinimizerDONLP2orthogonalproj.c

  using the minimization routine donlp2

  this is ugly as hell, but it works !


  (c) 2003-2006 Thomas Neff

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
#include "fmd/ConstraintJ2.h"
#include "fmd/Projection.h"
#include "numerics/donlp2.h"

#include "fmd/Hamiltonian.h"
#include "fmd/gradHamiltonian.h"

#include "MinimizerDONLP2orthogonalproj.h"

#include "numerics/fortranc.h"
#include "numerics/wignerd.h"
#include "misc/physics.h"
#include "misc/utils.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#endif

// boundary constraint for real part of a
#define MINIMUMA 0.0001
#define ALLOWEDJ2DEVIATION 0.1
#define ENERGYPENALTY 10000.0/hbc


// perform projection in body-fixed system
// to do that the body-fixed system is rotated into
// the lab system, the gradients have then to be
// rotated back

#define ORIENTED

// randomize the azimuthal angles to avoid
// falling into spurious minima

// #define RANDOMIZED


#ifdef ORIENTED
#include "fmd/SpatialOrientation.h"
#endif

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
  int j;
  int k;
  int par;
  angintegrationpara* projpar;
  int nangles;
  const Constraint* Const;
  int nconst;
  SlaterDet* Q0;
  double h0, n0;
  Para* q;
  Parameterization* P;
  SlaterDet* Q;
  SlaterDet* Qp;
  SlaterDet* Qpp;
  SlaterDetAux* X;
  gradSlaterDetAux* dX;
  gradSlaterDet *dH, *dN, *dh, *dn;
  gradSlaterDet* dHq0, *dNq0, *dHqq, *dNqq;
} Min;


static struct {
  SlaterDet* Qp;
  SlaterDet* Qpp;
  SlaterDetAux* X;
} Work;
  

#ifdef ORIENTED
static double lalpha0;
static double lbeta0;
static double lgamma0;

#define M_2PI (2*M_PI)
#define M_PI2 (M_PI*M_PI)

inline static double dabs(double x) { if (x<0) return -x; else return x; };

inline static double mod2pi(double a) { return (a - floor(a/M_2PI)*M_2PI); }

inline static double diffang(double a, double b)
{
  return (M_PI - dabs(M_PI-mod2pi(a-b)));
}

int orientationinverted(double alpha, double beta, double gamma,
			double lalpha, double lbeta, double lgamma)
{
  // if we are away from beta = pi/2 only look at beta
  // if we are close to beta = pi/2 also look at alpha
  if ((diffang(beta, M_PI_2) > 0.1 && 
       diffang(M_PI-beta, lbeta) < diffang(beta, lbeta)) ||
      (diffang(beta, M_PI_2) <= 0.1 && 
       diffang(M_PI-beta, lbeta)+diffang(alpha-M_PI, lalpha) < 
       diffang(beta, lbeta)+diffang(alpha, lalpha)))
    return 1;
  else
    return 0;
}


void invertangles(double* alpha, double* beta, double* gamma)
{
  *alpha = mod2pi(*alpha - M_PI);
  *beta  = M_PI - *beta;
  *gamma = mod2pi(M_PI - *gamma);
}
#endif


void initWork(const SlaterDet* Q)
{
  Work.Qp = malloc(sizeof(SlaterDet));
  Work.Qpp = malloc(sizeof(SlaterDet));
  Work.X = malloc(sizeof(SlaterDetAux));

  initSlaterDet(Q, Work.Qp);
  initSlaterDet(Q, Work.Qpp);
  initSlaterDetAux(Q, Work.X);
}


void calcprojectedHamiltonian(const SlaterDet* Q, 
			      const Interaction* Int,
			      int j, int k, int par, 
			      angintegrationpara* projpar,
			      double* hproj, double* nproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDet* Qpp = Work.Qpp;
  SlaterDetAux* X = Work.X;

  int n = projpar->n;

  double alpha, beta, gamma, w;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif
  complex double hme, nme, j2me;
  complex double H=0.0, N=0.0, J2=0.0;
   
  complex double hintr, nintr;

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0))
    invertangles(&alpha0, &beta0, &gamma0);
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);
#endif

  calcSlaterDetAuxod(Qp, Qp, X);

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  int i, ip;
  for (i=0; i<n; i++) {
    copySlaterDet(Qp, Qpp);

    getangintegrationpoint(i, projpar, &alpha, &beta, &gamma, &w);
#ifdef RANDOMIZED
    alpha += alpharan;
    gamma += gammaran;
#endif

    rotateSlaterDet(Qpp, alpha, beta, gamma);
    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      nme = X->ovlap;
      calcHamiltonianod(Int, Qp, Qpp, X, &hme);
      calcConstraintJ2od(Qp, Qpp, X, &j2me);

      H += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*hme;
      N += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*nme;
      J2 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*j2me;
    }
  }

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(N/nintr), 
	  creal(J2/N));

  *hproj = creal(H);
  *nproj = creal(N);
}


#ifdef MPI
void calcprojectedHamiltonianmpi(const SlaterDet* Q, 
				 const Interaction* Int,
				 int j, int k, int par, 
				 angintegrationpara* angpara,
				 double* hproj, double* nproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDetAux* X = Work.X;

  int nang = angpara->n;

#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double H=0.0, N=0.0, J2=0.0;
   
  complex double hintr, nintr;

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0)) {
    fprintf(stderr, "... angles inverted\n");
    invertangles(&alpha0, &beta0, &gamma0);
  }
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);
#endif

  calcSlaterDetAuxod(Qp, Qp, X);

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  int task = TASKHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Qp);
  BroadcastSlaterDet(Qp);

  int todo = nang;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[3];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;

  complex double h[2], n[2], j2[2];

  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getangintegrationpoint(next, angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weight[processor]); 

      for (int a=0; a<3; a++)
	projpar[a] = angle[processor][a];

      MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;

    } else {					// collect results from slaves

      MPI_Recv(h, 2, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGHAMILTONIANOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      MPI_Recv(n, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGOVLAPOD, MPI_COMM_WORLD, &status);
      MPI_Recv(j2, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGANGULARMOMENTUMOD, MPI_COMM_WORLD, &status);

      complex double w;
      w = weight[processor]*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,angle[processor][0],angle[processor][1],angle[processor][2]);

      H += w/2*(h[0] + par* h[1]);
      N += w/2*(n[0] + par* n[1]);
      J2 += w/2*(j2[0] + par* j2[1]);

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);
  }	

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(N/nintr), 
	  creal(J2/N));

  *hproj = creal(H);
  *nproj = creal(N);
}
#endif


#define SQR(x) ((x)*(x))

void calcorthprojectedHamiltonian(const SlaterDet* Q0, 
				  double H00, double N00,
				  const SlaterDet* Q, 
				  const Interaction* Int,
				  int j, int k, int par, 
				  angintegrationpara* projpar,
				  double* eintr, double* eproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDet* Qpp = Work.Qpp;
  SlaterDetAux* X = Work.X;

  int n = projpar->n;

  double alpha, beta, gamma, w;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif
  complex double hqqme, nqqme, j2qqme;
  complex double Hqq=0.0, Nqq=0.0, J2qq=0.0;

  complex double hq0me, nq0me, j2q0me;
  complex double Hq0=0.0, Nq0=0.0, J2q0=0.0;

  complex double H, N, J2;

  complex double hintr, nintr;

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0))
    invertangles(&alpha0, &beta0, &gamma0);
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);
#endif

  calcSlaterDetAuxod(Qp, Qp, X);

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  *eintr = creal(hintr/nintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  int i, ip;
  for (i=0; i<n; i++) {

    getangintegrationpoint(i, projpar, &alpha, &beta, &gamma, &w);
#ifdef RANDOMIZED
    alpha += alpharan;
    gamma += gammaran;
#endif

    copySlaterDet(Q0, Qpp);
    rotateSlaterDet(Qpp, alpha, beta, gamma);
    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      nq0me = X->ovlap;
      calcHamiltonianod(Int, Qp, Qpp, X, &hq0me);
      calcConstraintJ2od(Qp, Qpp, X, &j2q0me);

      Hq0 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*hq0me;
      Nq0 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*nq0me;
      J2q0 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*j2q0me;
    }

    copySlaterDet(Qp, Qpp);
    rotateSlaterDet(Qpp, alpha, beta, gamma);
    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      nqqme = X->ovlap;
      calcHamiltonianod(Int, Qp, Qpp, X, &hqqme);
      calcConstraintJ2od(Qp, Qpp, X, &j2qqme);

      Hqq += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*hqqme;
      Nqq += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*nqqme;
      J2qq += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*j2qqme;
    }
  }

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(Nqq/nintr), 
	  creal(J2qq/Nqq));

  H = Hqq - 2.0*creal(Hq0*conj(Nq0))/N00 + H00*SQR(cabs(Nq0)/N00);
  N = Nqq - SQR(cabs(Nq0))/N00;

  *eproj = creal(H/N);
}


#ifdef MPI
void calcorthprojectedHamiltonianmpi(const SlaterDet* Q0,
				     double H00, double N00,
				     const SlaterDet* Q, 
				     const Interaction* Int,
				     int j, int k, int par, 
				     angintegrationpara* angpara,
				     double* eintr, double* eproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDetAux* X = Work.X;

  int nang = angpara->n;

#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double Hqq=0.0, Nqq=0.0, J2qq=0.0;
  complex double Hq0=0.0, Nq0=0.0, J2q0=0.0;

  complex double H, N, J2;

  complex double hintr, nintr;

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0)) {
    fprintf(stderr, "... angles inverted\n");
    invertangles(&alpha0, &beta0, &gamma0);
  }
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);
#endif

  calcSlaterDetAuxod(Qp, Qp, X);

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  *eintr = creal(hintr/nintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  // code duplication for <Q| |Q0> and <Q| |Q> matrix elements

  {
  int task = TASKHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Qp);
  BroadcastSlaterDet(Q0);

  int todo = nang;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[3];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;

  complex double h[2], n[2], j2[2];

  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getangintegrationpoint(next, angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weight[processor]); 

      for (int a=0; a<3; a++)
	projpar[a] = angle[processor][a];

      MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;

    } else {					// collect results from slaves

      MPI_Recv(h, 2, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGHAMILTONIANOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      MPI_Recv(n, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGOVLAPOD, MPI_COMM_WORLD, &status);
      MPI_Recv(j2, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGANGULARMOMENTUMOD, MPI_COMM_WORLD, &status);

      complex double w;
      w = weight[processor]*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,angle[processor][0],angle[processor][1],angle[processor][2]);

      Hq0 += w/2*(h[0] + par* h[1]);
      Nq0 += w/2*(n[0] + par* n[1]);
      J2q0 += w/2*(j2[0] + par* j2[1]);

      working[processor] = 0;
      running--;
      done++;
    }

  }
   
  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);
  }	
  }

  {
  int task = TASKHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Qp);
  BroadcastSlaterDet(Qp);

  int todo = nang;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[3];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;

  complex double h[2], n[2], j2[2];

  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getangintegrationpoint(next, angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weight[processor]); 

      for (int a=0; a<3; a++)
	projpar[a] = angle[processor][a];

      MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;

    } else {					// collect results from slaves

      MPI_Recv(h, 2, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGHAMILTONIANOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      MPI_Recv(n, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGOVLAPOD, MPI_COMM_WORLD, &status);
      MPI_Recv(j2, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGANGULARMOMENTUMOD, MPI_COMM_WORLD, &status);

      complex double w;
      w = weight[processor]*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,angle[processor][0],angle[processor][1],angle[processor][2]);

      Hqq += w/2*(h[0] + par* h[1]);
      Nqq += w/2*(n[0] + par* n[1]);
      J2qq += w/2*(j2[0] + par* j2[1]);

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);
  }	
  }

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(Nqq/nintr), 
	  creal(J2qq/Nqq));

  H = Hqq - 2.0*creal(Hq0*conj(Nq0))/N00 + H00*SQR(cabs(Nq0)/N00);
  N = Nqq - SQR(cabs(Nq0))/N00;

  *eproj = creal(H/N);
}
#endif


void calcgradorthprojectedHamiltonian(const SlaterDet* Q, 
				      const Interaction* Int,
				      int j, int k, int par, 
				      angintegrationpara* projpar,
				      double* eintr, double* eproj, 
				      double* ovlap)
{
  SlaterDet* Qp = Min.Qp;
  SlaterDet* Qpp = Min.Qpp;
  SlaterDetAux* X = Min.X;
  gradSlaterDetAux* dX = Min.dX;
  gradSlaterDet* dH = Min.dH;
  gradSlaterDet* dN = Min.dN;
  gradSlaterDet* dh = Min.dh;
  gradSlaterDet* dn = Min.dn;

  SlaterDet* Q0 = Min.Q0;
  double H00 = Min.h0; double N00 = Min.n0;
  gradSlaterDet *dHqq = Min.dHqq;
  gradSlaterDet *dNqq = Min.dNqq;
  gradSlaterDet *dHq0 = Min.dHq0;
  gradSlaterDet *dNq0 = Min.dNq0;

  int n = projpar->n;
  double alpha, beta, gamma, w;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double Hq0=0.0, Nq0=0.0;
  complex double Hqq=0.0, Nqq=0.0;
  complex double H=0.0, N=0.0;
  zerogradSlaterDet(dH);
  zerogradSlaterDet(dN);
  zerogradSlaterDet(dHqq);
  zerogradSlaterDet(dNqq);
  zerogradSlaterDet(dHq0);
  zerogradSlaterDet(dNq0);

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0)) {
    fprintf(stderr, "... angles inverted\n");
    invertangles(&alpha0, &beta0, &gamma0);
  }
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);

  lalpha0 = alpha0; lbeta0 = beta0; lgamma0 = gamma0;
#endif

  complex double hintr, nintr;
  calcSlaterDetAuxod(Qp, Qp, X);
  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  *eintr= creal(hintr/nintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  int i, ip;
  for (i=0; i<n; i++) {

    getangintegrationpoint(i, projpar, &alpha, &beta, &gamma, &w);
#ifdef RANDOMIZED
    alpha += alpharan;
    gamma += gammaran;
#endif

    copySlaterDet(Q0, Qpp);
    rotateSlaterDet(Qpp, alpha, beta, gamma);

    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      calcgradSlaterDetAuxod(Qp, Qpp, X, dX);
      calcgradOvlapod(Qp, Qpp, X, dX, dn);
      dn->val = X->ovlap;
      calcgradHamiltonianod(Int, Qp, Qpp, X, dX, dh);
 
      Hq0 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*dh->val;
      Nq0 += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*dn->val;

      addmulttogradSlaterDet(dHq0, dh, w/2*(ip ? par : +1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma));
      addmulttogradSlaterDet(dNq0, dn, w/2*(ip ? par : +1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma));
    }

    copySlaterDet(Qp, Qpp);
    rotateSlaterDet(Qpp, alpha, beta, gamma);

    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      calcgradSlaterDetAuxod(Qp, Qpp, X, dX);
      calcgradOvlapod(Qp, Qpp, X, dX, dn);
      dn->val = X->ovlap;
      calcgradHamiltonianod(Int, Qp, Qpp, X, dX, dh);
 
      Hqq += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*dh->val;
      Nqq += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma)*dn->val;

      addmulttogradSlaterDet(dHqq, dh, w/2*(ip ? par : +1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma));
      addmulttogradSlaterDet(dNqq, dn, w/2*(ip ? par : +1)*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,alpha,beta,gamma));
    }
  }

  H = Hqq - 2.0*creal(Hq0*conj(Nq0))/N00 + H00*SQR(cabs(Nq0)/N00);
  N = Nqq - SQR(cabs(Nq0))/N00;

  addmulttogradSlaterDet(dH, dHqq, 1.0);
  addmulttogradSlaterDet(dN, dNqq, 1.0);

  addmulttogradSlaterDet(dH, dHq0, -conj(Nq0)/N00);
  addmulttogradSlaterDet(dH, dNq0, H00*conj(Nq0)/SQR(N00)- conj(Hq0)/N00);
  addmulttogradSlaterDet(dN, dNq0, -conj(Nq0)/N00);

  *ovlap = cabs(Nq0)/sqrt(creal(Nqq*N00));

  *eproj = creal(H/N);
  multgradSlaterDet(dH, 1.0/N);
  addmulttogradSlaterDet(dH, dN, -H/(N*N));

#ifdef ORIENTED
  rotategradSlaterDet(dH, alpha0, beta0, gamma0);
#endif
}				  


#ifdef MPI
void calcgradorthprojectedHamiltonianmpi(const SlaterDet* Q, 
					 const Interaction* Int,
					 int j, int k, int par, 
					 angintegrationpara* angpara,
					 double* eintr, double* eproj,
					 double* ovlap)
{
  SlaterDet* Qp = Min.Qp;
  SlaterDetAux* X = Min.X;
  gradSlaterDet* dH = Min.dH;
  gradSlaterDet* dN = Min.dN;
  gradSlaterDet* dh = Min.dh;
  gradSlaterDet* dn = Min.dn;

  SlaterDet* Q0 = Min.Q0;
  double H00 = Min.h0; double N00 = Min.n0;
  gradSlaterDet *dHqq = Min.dHqq;
  gradSlaterDet *dNqq = Min.dNqq;
  gradSlaterDet *dHq0 = Min.dHq0;
  gradSlaterDet *dNq0 = Min.dNq0;

  int nang = angpara->n;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double Hq0=0.0, Nq0=0.0;
  complex double Hqq=0.0, Nqq=0.0;
  complex double H=0.0, N=0.0;
  zerogradSlaterDet(dH);
  zerogradSlaterDet(dN);
  zerogradSlaterDet(dHqq);
  zerogradSlaterDet(dNqq);
  zerogradSlaterDet(dHq0);
  zerogradSlaterDet(dNq0);

  copySlaterDet(Q, Qp);
#ifdef ORIENTED
  calcSlaterDetAuxod(Qp, Qp, X);
  calcOrientedOrientation(Qp, X, &alpha0, &beta0, &gamma0);
  if (orientationinverted(alpha0, beta0, gamma0, lalpha0, lbeta0, lgamma0)) {
    fprintf(stderr, "... angles inverted\n");
    invertangles(&alpha0, &beta0, &gamma0);
  }
  fprintf(stderr, "... (alpha, beta, gamma) = (%5.2f, %5.2f, %5.2f)\n",
	  alpha0*180/M_PI, beta0*180/M_PI, gamma0*180/M_PI); 
  rotateSlaterDet(Qp, -gamma0, -beta0, -alpha0);

  lalpha0 = alpha0; lbeta0 = beta0; lgamma0 = gamma0;
#endif

  complex double hintr, nintr;
  calcSlaterDetAuxod(Qp, Qp, X);
  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  *eintr= creal(hintr/nintr);

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  {
  int task = TASKGRADHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Qp);
  BroadcastSlaterDet(Q0);

  int todo = nang;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[3];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;


  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getangintegrationpoint(next, angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weight[processor]); 
      for (int k=0; k<3; k++)
	projpar[k] = angle[processor][k];

      MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;

    } else {					// collect results from slaves

      // wait for first value, identify slave, then recv all the other stuff from this slave

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);

      processor = status.MPI_SOURCE;
      complex double w;
      w = weight[processor]*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,angle[processor][0],angle[processor][1],angle[processor][2]);

      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

      Hq0 += w/2* dh->val;
      addmulttogradSlaterDet(dHq0, dh, 0.5*w);    

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

      Hq0 += w/2*par* dh->val;
      addmulttogradSlaterDet(dHq0, dh, 0.5*par*w);

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);
    
      Nq0 += w/2* dn->val;
      addmulttogradSlaterDet(dNq0, dn, 0.5*w);    

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);

      Nq0 += w/2*par* dn->val;
      addmulttogradSlaterDet(dNq0, dn, 0.5*par*w);    

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);
  }	

  }

  {
  int task = TASKGRADHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Qp);
  BroadcastSlaterDet(Qp);

  int todo = nang;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[3];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;


  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getangintegrationpoint(next, angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weight[processor]); 
      for (int k=0; k<3; k++)
	projpar[k] = angle[processor][k];

      MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;

    } else {					// collect results from slaves

      // wait for first value, identify slave, then recv all the other stuff from this slave

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);

      processor = status.MPI_SOURCE;
      complex double w;
      w = weight[processor]*(j+1)/(8*M_PI2)*Djmkstar(j,k,k,angle[processor][0],angle[processor][1],angle[processor][2]);

      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

      Hqq += w/2* dh->val;
      addmulttogradSlaterDet(dHqq, dh, 0.5*w);    

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

      Hqq += w/2*par* dh->val;
      addmulttogradSlaterDet(dHqq, dh, 0.5*par*w);

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);
    
      Nqq += w/2* dn->val;
      addmulttogradSlaterDet(dNqq, dn, 0.5*w);    

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);

      Nqq += w/2*par* dn->val;
      addmulttogradSlaterDet(dNqq, dn, 0.5*par*w);    

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 3, MPI_DOUBLE, processor, TAGPROJECT3, MPI_COMM_WORLD);
  }	

  }

  H = Hqq - 2.0*creal(Hq0*conj(Nq0))/N00 + H00*SQR(cabs(Nq0)/N00);
  N = Nqq - SQR(cabs(Nq0))/N00;

  addmulttogradSlaterDet(dH, dHqq, 1.0);
  addmulttogradSlaterDet(dN, dNqq, 1.0);

  addmulttogradSlaterDet(dH, dHq0, -conj(Nq0)/N00);
  addmulttogradSlaterDet(dH, dNq0, H00*conj(Nq0)/SQR(N00)- conj(Hq0)/N00);
  addmulttogradSlaterDet(dN, dNq0, -conj(Nq0)/N00);

  *ovlap = cabs(Nq0)/sqrt(creal(Nqq*N00));

  *eproj = creal(H/N);
  multgradSlaterDet(dH, 1.0/N);
  addmulttogradSlaterDet(dH, dN, -H/(N*N));

#ifdef ORIENTED
  rotategradSlaterDet(dH, alpha0, beta0, gamma0);
#endif
}				  
#endif



static void copyxtopara(const double* x, Para* q)
{
  int i;
  for (i=0; i<q->n; i++)
    q->x[i] = x[i];
}


void MinimizeDONLP2orthproj(const Interaction* Int, int j, int k, int par,
			    angintegrationpara* projpar,
			    SlaterDet* Q0, double h0, double n0,
			    const Constraint* Const, int nconst,
			    Parameterization* P, Para* q,
			    int maxsteps, int log, const char* logfile)
{
#ifndef MPI
  // handler for INT and TERM signals
  signal(SIGINT, catchterminate);
  signal(SIGTERM, catchterminate);
#endif

  Min.maxsteps = maxsteps;
  Min.log = log;
  Min.logfile = logfile;
  Min.Int = Int;
  Min.j = j;
  Min.k = k;
  Min.par = par;
  Min.projpar = projpar;
  Min.Const = Const;
  Min.nconst = nconst;
  Min.Q0 = Q0;
  Min.h0 = h0;
  Min.n0 = n0;
  Min.q = q;
  Min.P = P;
  Min.Q = malloc(sizeof(SlaterDet));
  Min.Qp = malloc(sizeof(SlaterDet));
  Min.Qpp = malloc(sizeof(SlaterDet));
  Min.X = malloc(sizeof(SlaterDetAux));
  Min.dX = malloc(sizeof(gradSlaterDetAux));
  Min.dH = malloc(sizeof(gradSlaterDet));
  Min.dN = malloc(sizeof(gradSlaterDet));
  Min.dh = malloc(sizeof(gradSlaterDet));
  Min.dn = malloc(sizeof(gradSlaterDet));
  Min.dHq0 = malloc(sizeof(gradSlaterDet));
  Min.dNq0 = malloc(sizeof(gradSlaterDet));
  Min.dHqq = malloc(sizeof(gradSlaterDet));
  Min.dNqq = malloc(sizeof(gradSlaterDet));
  Min.P->ParainitSlaterDet(q, Min.Q);
  Min.P->ParainitSlaterDet(q, Min.Qp);
  Min.P->ParainitSlaterDet(q, Min.Qpp);
  initSlaterDetAux(Min.Q, Min.X);
  initgradSlaterDetAux(Min.Q, Min.dX);
  initgradSlaterDet(Min.Q, Min.dH);
  initgradSlaterDet(Min.Q, Min.dN);
  initgradSlaterDet(Min.Q, Min.dh);
  initgradSlaterDet(Min.Q, Min.dn);
  initgradSlaterDet(Min.Q, Min.dHq0);
  initgradSlaterDet(Min.Q, Min.dNq0);
  initgradSlaterDet(Min.Q, Min.dHqq);
  initgradSlaterDet(Min.Q, Min.dNqq);
  

#ifndef MPI
  // regular program execution or catched signal ?
  if (!setjmp(env))
#endif
    FORTRAN(donlp2)();

  fprintf(stderr, "donlp2 terminated because of criterium %d\n",
	  (int) FORTRAN(o8itin).optite + 11);

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

  double eintr, eproj;
#ifdef MPI
  calcorthprojectedHamiltonianmpi(Min.Q0, Min.h0, Min.n0,
				  Min.Q,
				  Min.Int,
				  Min.j, Min.k, Min.par, Min.projpar,
				  &eintr, &eproj);
#else
  calcorthprojectedHamiltonian(Min.Q0, Min.h0, Min.n0,
			       Min.Q,
			       Min.Int,
			       Min.j, Min.k, Min.par, Min.projpar,
			       &eintr, &eproj);
#endif

  *fx = eproj;
  FORTRAN(o8cnt).icf++;

  fprintf(stderr, "\t\t\tme      E = %8.3f MeV,   Eintr = %8.3f MeV\n", 
	  hbc*eproj, hbc*eintr);

}


void FORTRAN(egradf)(const double* x, double* gradf)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif

  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);

  double eintr, eproj, ovlap;

#ifdef MPI
  calcgradorthprojectedHamiltonianmpi(Min.Q,
				      Min.Int,
				      Min.j, Min.k, Min.par, Min.projpar,
				      &eintr, &eproj, &ovlap);
#else
  calcgradorthprojectedHamiltonian(Min.Q,
				   Min.Int,
				   Min.j, Min.k, Min.par, Min.projpar,
				   &eintr, &eproj, &ovlap);
#endif

  Min.P->ParaprojectgradSlaterDet(Min.q, Min.dH, gradf);
  FORTRAN(o8cnt).icgf++;
  fprintf(stderr, "grad %3d: \tE = %8.3f MeV,   Eintr = %8.3f MeV,  ovlap = %5.3f\n", 
	  FORTRAN(o8cnt).icgf, hbc*eproj, hbc*eintr, ovlap);


  if (Min.log && !(FORTRAN(o8cnt).icgf % Min.log)) {
    char logfname[255];
    sprintf(logfname, "%s.minproj.%03d.log", Min.logfile, FORTRAN(o8cnt).icgf);
    FILE* logfp;
    if (!(logfp = fopen(logfname, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", logfname);
    } else {
      fprintinfo(logfp);
      fprintf(logfp, "# step %3d: \tE = %8.3f MeV,   Eintr = %8.3f MeV,  ovlap = %5.3f\n", 
	      FORTRAN(o8cnt).icgf, hbc*eproj, hbc*eintr, ovlap);
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
    *gxi = x[(*i-1)*NGAUSS+APOS]-MINIMUMA;
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
  FORTRAN(o8par).tau0 = 20.0;
  FORTRAN(o8par).del0 = 1.0;
  FORTRAN(o8rst).nreset = 20;
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
