/**

  \file MinimizerDONLP2vap.c

  using the minimization routine donlp2

  this is ugly as hell, but it works !


  (c) 2003-2007 Thomas Neff

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

#include "MinimizerDONLP2vap.h"

#include "numerics/fortranc.h"
#include "numerics/wignerd.h"
#include "numerics/lapack.h"
#include "numerics/cmat.h"
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

// K-mixing threshold
#define THRESHKMIX 0.01


// perform projection in body-fixed system
// to do that the body-fixed system is rotated into
// the lab system, the gradients have then to be
// rotated back

#define ORIENTED


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
  int par;
  angintegrationpara* projpar;
  int nangles;
  const Constraint* Const;
  int nconst;
  Para* q;
  Parameterization* P;
  SlaterDet* Q;
  SlaterDet* Qp;
  SlaterDet* Qpp;
  SlaterDetAux* X;
  gradSlaterDetAux* dX;
  gradSlaterDet *dh, *dn, *dH, *dN, *dhproj, *dnproj;
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


// solve generalized eigenvalue problem, 
// return lowest eigenvalue and corresponding eigenvector
static int geneigproblem(int n, complex double* A, complex double* B, double thresh, 
			 complex double* a, complex double* c)
{
  const char jobu='O', jobvt='A';

  complex double U[n*n];
  complex double Vh[n*n];
  double S[n];

  const int lwork=5*n;
  complex double work[lwork];
  double rwork[lwork];

  int info;
  
  // perform SVD of overlap matrix
  copycmat(n, B, U);
  FORTRAN(zgesvd)(&jobu, &jobvt, &n, &n, U, &n, S, 
		  NULL, &n, Vh, &n, work, &lwork, rwork, &info);

  // for zero matrices
  if (S[0] == 0.0) {
    return -1;
  }

  // dimension of subspace
  int m=0;

  while (m<n && S[m]/S[0] > thresh)
    m++;

  // if our subspace is empty nothing to do any more

  if (!m) {
    return -1;
  }

  // project into subspace
  complex double AV[m*n];
  complex double invNA[m*m];
  int k,l,i,j;

  for (i=0; i<n; i++)
    for (l=0; l<m; l++) {
      AV[i+l*n] = 0.0;
      for (j=0; j<n; j++)
	AV[i+l*n] += A[i+j*n]*conj(Vh[l+j*n]);
    }
	
  for (l=0; l<m; l++)
    for (k=0; k<m; k++) {
      invNA[k+l*m] = 0.0;
      for (i=0; i<n; i++)
	invNA[k+l*m] += 1/S[k]*conj(U[i+k*n])*AV[i+l*n];
      }	

  // solve eigenvalue problem in subspace

  const char jobvl='N', jobvr='V';
  complex double vb[m];
  complex double Vb[m*m];

  FORTRAN(zgeev)(&jobvl, &jobvr, &m, invNA, &m, vb, NULL, &m, Vb, &m,
		 work, &lwork, rwork, &info);

  // find smalles eigenvalue
  complex double vmin = vb[0];
  int imin = 0;

  for (i=1; i<m; i++)
    if (creal(vb[i]) < creal(vmin)) {
      vmin = vb[i];
      imin = i;
    }

  // embed solution in full space

  for (i=0; i<n; i++) {
    c[i] = 0.0;
    for (l=0; l<m; l++)
      c[i] += conj(Vh[l+i*n])*Vb[l+imin*m];
    }

  return 0;
}


static complex double expect(int n, complex double* A, complex double* c)
{
  int i,j;
  complex double res = 0.0;

  for (j=0; j<n; j++)
    for (i=0; i<n; i++)
      res += conj(c[i])*A[i+j*n]*c[j];

  return res;
}


static void gradexpect(int n, gradSlaterDet* dA, complex double* c,
		       gradSlaterDet* da)
{
  int i,j;

  zerogradSlaterDet(da);

  for (j=0; j<n; j++)
    for (i=0; i<n; i++)
      addmulttogradSlaterDet(da, &dA[i+j*n], conj(c[i])*c[j]);
}

void calcprojectedHamiltonian(const SlaterDet* Q, 
			      const Interaction* Int,
			      int j, int par, 
			      angintegrationpara* angpara,
			      double* eintr, double* eproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDet* Qpp = Work.Qpp;
  SlaterDetAux* X = Work.X;

  int n = angpara->n;

  double alpha, beta, gamma, w;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif
  complex double hme, nme, j2me;
  complex double H[(j+1)*(j+1)], N[(j+1)*(j+1)], J2[(j+1)*(j+1)];

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      J2[idxjmk(j,k,kp)] = 0.0;
    }
   
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

  reinitangintegration(Qp, Qp, 0, 0, angpara);
  calcSlaterDetAuxod(Qp, Qp, X);

  if (isnan(X->ovlap)) {
    *eintr = ENERGYPENALTY;
    *eproj = ENERGYPENALTY;
    return;
  }

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  *eintr = hintr/nintr;

#ifdef RANDOMIZED
  double alpharan, gammaran;
  alpharan = 0.005*2*M_PI*rand()/RAND_MAX;
  gammaran = 0.005*2*M_PI*rand()/RAND_MAX;
#endif

  int i, ip;
  for (i=0; i<n; i++) {
    copySlaterDet(Qp, Qpp);

    getangintegrationpoint(i, angpara, &alpha, &beta, &gamma, &w);
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

      for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  H[idxjmk(j,k,kp)] += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,alpha,beta,gamma)*hme;
	  N[idxjmk(j,k,kp)] += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,alpha,beta,gamma)*nme;
	  J2[idxjmk(j,k,kp)] += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,alpha,beta,gamma)*j2me;
	}
    }
  }

  complex double c[j+1];
  complex double e;
  geneigproblem(j+1, H, N, THRESHKMIX, &e, c);

  complex double hproj, nproj, j2proj;
  hproj = expect(j+1, H, c);
  nproj = expect(j+1, N, c);
  j2proj = expect(j+1, J2, c);

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(nproj/nintr), 
	  creal(j2proj/nproj));

  if (dabs(j2proj/nproj - 0.5*j*(0.5*j+1.0)) > ALLOWEDJ2DEVIATION) {
    *eintr = ENERGYPENALTY;
    *eproj = ENERGYPENALTY;
    return;
  }

  *eproj = creal(hproj/nproj);
}


#ifdef MPI
void calcprojectedHamiltonianmpi(const SlaterDet* Q, 
				 const Interaction* Int,
				 int j, int par, 
				 angintegrationpara* angpara,
				 double* eintr, double* eproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDetAux* X = Work.X;

  int nang = angpara->n;

#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double H[(j+1)*(j+1)], N[(j+1)*(j+1)], J2[(j+1)*(j+1)];

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      J2[idxjmk(j,k,kp)] = 0.0;
    }
   
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

  reinitangintegration(Qp, Qp, 0, 0, angpara);
  calcSlaterDetAuxod(Qp, Qp, X);

  if (isnan(X->ovlap)) {
    *eintr = ENERGYPENALTY;
    *eproj = ENERGYPENALTY;
    return;
  }

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);
  *eintr = hintr/nintr;

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

  complex double hme[2], nme[2], j2me[2];

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

      MPI_Recv(hme, 2, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGHAMILTONIANOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      MPI_Recv(nme, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGOVLAPOD, MPI_COMM_WORLD, &status);
      MPI_Recv(j2me, 2, MPI_DOUBLE_COMPLEX, 
	       processor, TAGANGULARMOMENTUMOD, MPI_COMM_WORLD, &status);

      complex double w;
      for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  w = weight[processor]*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,angle[processor][0],angle[processor][1],angle[processor][2]);

	  H[idxjmk(j,k,kp)] += w/2*(hme[0] + par* hme[1]);
	  N[idxjmk(j,k,kp)] += w/2*(nme[0] + par* nme[1]);
	  J2[idxjmk(j,k,kp)] += w/2*(j2me[0] + par* j2me[1]);
	}

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

  complex double c[j+1];
  complex double e;
  geneigproblem(j+1, H, N, THRESHKMIX, &e, c);

  complex double hproj, nproj, j2proj;
  hproj = expect(j+1, H, c);
  nproj = expect(j+1, N, c);
  j2proj = expect(j+1, J2, c);

  fprintf(stderr, "... n: %9.7f, J2: %6.3f\n", 
	  creal(nproj/nintr), 
	  creal(j2proj/nproj));

  if (dabs(j2proj/nproj - 0.5*j*(0.5*j+1.0)) > ALLOWEDJ2DEVIATION) {
    *eintr = ENERGYPENALTY;
    *eproj = ENERGYPENALTY;
    return;
  }

  *eproj = creal(hproj/nproj);
}
#endif


void calcgradprojectedHamiltonian(const SlaterDet* Q, 
				  const Interaction* Int,
				  int j, int par, 
				  angintegrationpara* angpara,
				  double* eintr, double* eproj)
{
  SlaterDet* Qp = Min.Qp;
  SlaterDet* Qpp = Min.Qpp;
  SlaterDetAux* X = Min.X;
  gradSlaterDetAux* dX = Min.dX;
  gradSlaterDet* dh = Min.dh;
  gradSlaterDet* dn = Min.dn;
  gradSlaterDet* dH = Min.dH;
  gradSlaterDet* dN = Min.dN;
  gradSlaterDet* dhproj = Min.dhproj;
  gradSlaterDet* dnproj = Min.dnproj;

  int n = angpara->n;
  double alpha, beta, gamma, w;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double H[(j+1)*(j+1)], N[(j+1)*(j+1)];

  int k, kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; N[idxjmk(j,k,kp)] = 0.0; 
      zerogradSlaterDet(&dH[idxjmk(j,k,kp)]);
      zerogradSlaterDet(&dN[idxjmk(j,k,kp)]);
    }

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

  reinitangintegration(Q, Q, 0, 0, angpara);

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
    copySlaterDet(Qp, Qpp);

    getangintegrationpoint(i, angpara, &alpha, &beta, &gamma, &w);
#ifdef RANDOMIZED
    alpha += alpharan;
    gamma += gammaran;
#endif

    rotateSlaterDet(Qpp, alpha, beta, gamma);

    for (ip=0; ip<=1; ip++) {

      if (ip) invertSlaterDet(Qpp);
      calcSlaterDetAuxod(Qp, Qpp, X);
      calcgradSlaterDetAuxod(Qp, Qpp, X, dX);
      calcgradOvlapod(Qp, Qpp, X, dX, dn);
      dn->val = X->ovlap;
      calcgradHamiltonianod(Int, Qp, Qpp, X, dX, dh);
 
      for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  H[idxjmk(j,k,kp)] += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,alpha,beta,gamma)*dh->val;
	  N[idxjmk(j,k,kp)] += w/2*(ip ? par : 1)*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,alpha,beta,gamma)*dn->val;

	  addmulttogradSlaterDet(&dH[idxjmk(j,k,kp)], dh, w/2*(ip ? par : +1)*
				 (j+1)/(8*M_PI2)*Djmkstar(j,k,kp,alpha,beta,gamma));
	  addmulttogradSlaterDet(&dN[idxjmk(j,k,kp)], dn, w/2*(ip ? par : +1)*
				 (j+1)/(8*M_PI2)*Djmkstar(j,k,kp,alpha,beta,gamma));
	}
    }
  }

  complex double c[j+1];
  complex double e;
  geneigproblem(j+1, H, N, THRESHKMIX, &e, c);

  complex double hproj, nproj;
  hproj = expect(j+1, H, c);
  nproj = expect(j+1, N, c);

  *eproj = creal(hproj/nproj);

  gradexpect(j+1, dH, c, dhproj);
  gradexpect(j+1, dN, c, dnproj);

  multgradSlaterDet(dhproj, 1.0/nproj);
  addmulttogradSlaterDet(dhproj, dnproj, -hproj/(nproj*nproj));

#ifdef ORIENTED
  rotategradSlaterDet(dhproj, alpha0, beta0, gamma0);
#endif
}				  


#ifdef MPI
void calcgradprojectedHamiltonianmpi(const SlaterDet* Q, 
				     const Interaction* Int,
				     int j, int par, 
				     angintegrationpara* angpara,
				     double* eintr, double* eproj)
{
  SlaterDet* Qp = Min.Qp;
  SlaterDetAux* X = Min.X;
  gradSlaterDet* dh = Min.dh;
  gradSlaterDet* dn = Min.dn;
  gradSlaterDet* dH = Min.dH;
  gradSlaterDet* dN = Min.dN;
  gradSlaterDet* dhproj = Min.dhproj;
  gradSlaterDet* dnproj = Min.dnproj;

  int nang = angpara->n;
#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

  complex double H[(j+1)*(j+1)], N[(j+1)*(j+1)];

  int k, kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; N[idxjmk(j,k,kp)] = 0.0; 
      zerogradSlaterDet(&dH[idxjmk(j,k,kp)]);
      zerogradSlaterDet(&dN[idxjmk(j,k,kp)]);
    }

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

  reinitangintegration(Q, Q, 0, 0, angpara);

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

      complex double w;

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);

      processor = status.MPI_SOURCE;

      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

     for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  w = weight[processor]*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,angle[processor][0],angle[processor][1],angle[processor][2]);

	  H[idxjmk(j,k,kp)] += w/2* dh->val;
	  addmulttogradSlaterDet(&dH[idxjmk(j,k,kp)], dh, w/2);
	}    

      MPI_Recv(&dh->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dh->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD, &status);

     for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  w = weight[processor]*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,angle[processor][0],angle[processor][1],angle[processor][2]);

	  H[idxjmk(j,k,kp)] += w/2*par* dh->val;
	  addmulttogradSlaterDet(&dH[idxjmk(j,k,kp)], dh, w/2*par);
	}    

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);

     for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  w = weight[processor]*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,angle[processor][0],angle[processor][1],angle[processor][2]);

	  N[idxjmk(j,k,kp)] += w/2* dn->val;
	  addmulttogradSlaterDet(&dN[idxjmk(j,k,kp)], dn, w/2);
	}        

      MPI_Recv(&dn->val, 1, MPI_DOUBLE_COMPLEX, processor, TAGGRADOVLAPODVAL, MPI_COMM_WORLD, &status);
      MPI_Recv(dn->gradval, Q->ngauss*sizeof(gradGaussian),
	       MPI_BYTE, processor, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD, &status);

     for (kp=-j; kp<=j; kp=kp+2)
	for (k=-j; k<=j; k=k+2) {
	  w = weight[processor]*(j+1)/(8*M_PI2)*
	    Djmkstar(j,k,kp,angle[processor][0],angle[processor][1],angle[processor][2]);

	  N[idxjmk(j,k,kp)] += w/2*par* dn->val;
	  addmulttogradSlaterDet(&dN[idxjmk(j,k,kp)], dn, w/2*par);
	}        

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

  complex double c[j+1];
  complex double e;
  geneigproblem(j+1, H, N, THRESHKMIX, &e, c);

  complex double hproj, nproj;
  hproj = expect(j+1, H, c);
  nproj = expect(j+1, N, c);

  *eproj = creal(hproj/nproj);

  gradexpect(j+1, dH, c, dhproj);
  gradexpect(j+1, dN, c, dnproj);

  multgradSlaterDet(dhproj, 1.0/nproj);
  addmulttogradSlaterDet(dhproj, dnproj, -hproj/(nproj*nproj));

#ifdef ORIENTED
  rotategradSlaterDet(dhproj, alpha0, beta0, gamma0);
#endif
}				  
#endif



static void copyxtopara(const double* x, Para* q)
{
  int i;
  for (i=0; i<q->n; i++)
    q->x[i] = x[i];
}


void MinimizeDONLP2vap(const Interaction* Int, int j, int par,
		       angintegrationpara* projpar,
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
  Min.par = par;
  Min.projpar = projpar;
  Min.Const = Const;
  Min.nconst = nconst;
  Min.q = q;
  Min.P = P;
  Min.Q = malloc(sizeof(SlaterDet));
  Min.Qp = malloc(sizeof(SlaterDet));
  Min.Qpp = malloc(sizeof(SlaterDet));
  Min.X = malloc(sizeof(SlaterDetAux));
  Min.dX = malloc(sizeof(gradSlaterDetAux));
  Min.dh = malloc(sizeof(gradSlaterDet));
  Min.dn = malloc(sizeof(gradSlaterDet));
  Min.dH = malloc((j+1)*(j+1)*sizeof(gradSlaterDet));
  Min.dN = malloc((j+1)*(j+1)*sizeof(gradSlaterDet));
  Min.dhproj = malloc(sizeof(gradSlaterDet));
  Min.dnproj = malloc(sizeof(gradSlaterDet));
  Min.P->ParainitSlaterDet(q, Min.Q);
  Min.P->ParainitSlaterDet(q, Min.Qp);
  Min.P->ParainitSlaterDet(q, Min.Qpp);
  initSlaterDetAux(Min.Q, Min.X);
  initgradSlaterDetAux(Min.Q, Min.dX);
  initgradSlaterDet(Min.Q, Min.dh);
  initgradSlaterDet(Min.Q, Min.dn);
  for (int i=0; i<(j+1)*(j+1); i++) {
    initgradSlaterDet(Min.Q, &Min.dH[i]);
    initgradSlaterDet(Min.Q, &Min.dN[i]);
  }
  initgradSlaterDet(Min.Q, Min.dhproj);
  initgradSlaterDet(Min.Q, Min.dnproj);

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
  calcprojectedHamiltonianmpi(Min.Q,
			   Min.Int,
			   Min.j, Min.par, Min.projpar,
			   &eintr, &eproj);
#else
  calcprojectedHamiltonian(Min.Q,
			   Min.Int,
			   Min.j, Min.par, Min.projpar,
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

  double eintr, eproj;

#ifdef MPI
  calcgradprojectedHamiltonianmpi(Min.Q,
			       Min.Int,
			       Min.j, Min.par, Min.projpar,
			       &eintr, &eproj);
#else
  calcgradprojectedHamiltonian(Min.Q,
			       Min.Int,
			       Min.j, Min.par, Min.projpar,
			       &eintr, &eproj);
#endif

  Min.P->ParaprojectgradSlaterDet(Min.q, Min.dhproj, gradf);
  FORTRAN(o8cnt).icgf++;
  fprintf(stderr, "grad %3d: \tE = %8.3f MeV,   Eintr = %8.3f MeV\n", 
	  FORTRAN(o8cnt).icgf, hbc*eproj, hbc*eintr);


  if (Min.log && !(FORTRAN(o8cnt).icgf % Min.log)) {
    char logfname[255];
    sprintf(logfname, "%s.minvap.%03d.log", Min.logfile, FORTRAN(o8cnt).icgf);
    FILE* logfp;
    if (!(logfp = fopen(logfname, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", logfname);
    } else {
      fprintinfo(logfp);
      fprintf(logfp, "# step %3d: \tE = %8.3f MeV,   Eintr = %8.3f MeV\n", 
	      FORTRAN(o8cnt).icgf, hbc*eproj, hbc*eintr);
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
