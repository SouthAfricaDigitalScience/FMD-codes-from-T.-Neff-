/**

  \file MinimizerDONLP2multivapp.c

  Minimize Eproj + alpha Eintr in order to improve overlap of projected and intrinsic state

  using the minimization routine donlp2

  this is ugly as hell, but it works !


  (c) 2003-2008 Thomas Neff

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

#include "MinimizerDONLP2vapp.h"

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
  int ival;
  double minnormkmix;
  double threshkmix;
  double alpha;
  angintegrationpara* projpar;
  int nangles;
  const Constraint* Const;
  int nconst;
  Para* q;
  Parameterization* P;
  SlaterDet* Q;
} Min;


static struct {
  complex double* V  ;       // K-mixing eigenstates
  int dim;                   // dimension
  complex double* H;         // projected Hamiltonian matrix elements for 
  complex double* N;         // projected Overlap matrix elements
  complex double* J2;

  gradSlaterDet* dH;
  gradSlaterDet* dN;

  SlaterDet* Qp;
  SlaterDet* Qpp;
  SlaterDetAux* X;
  gradSlaterDetAux* dX;

  gradSlaterDet* dhproj;
  gradSlaterDet* dnproj;
} Work;
  

#define M_2PI (2*M_PI)
#define M_PI2 (M_PI*M_PI)

inline static double dabs(double x) { if (x<0) return -x; else return x; };

#ifdef ORIENTED
static double lalpha0;
static double lbeta0;
static double lgamma0;

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


// perform K-mixing calculation and sort out "good" states

static int cmpmerit(void* ap, void* bp) 
{
  double a=*(double*) ap;
  double b=*(double*) bp;
  
  return (a<=b ? (a<b ? 1 : 0) : -1);
}


void kmixingstates(complex double* H, complex double* N, int j, double thresh, 
                   double* egood, complex double* Vgood, int* ngood)
{
  complex double e[j+1];
  complex double V[(j+1)*(j+1)];
  int d;

  generalizedeigensystem(H, N, j+1, thresh, e, V, &d);

  // we got d eigenstates, filter out the noise

  // sort eigenstates according to energy and norm

  // calculate norms
  complex double N2[(j+1)*(j+1)];
  multcmat(N, N, N2, j+1);

  double n, n2, norm[d];
  int i, k, l;

  for (i=0; i<d; i++) {
    n=0.0; n2=0.0;
    for (l=0; l<j+1; l++)
      for (k=0; k<j+1; k++) {
        n += conj(V[k+i*(j+1)])*N[k+l*(j+1)]*V[l+i*(j+1)];
        n2 += conj(V[k+i*(j+1)])*N2[k+l*(j+1)]*V[l+i*(j+1)];
      }
    norm[i] = n2/n;
  }

  double maxnorm = 0.0;
  for (i=0; i<d; i++)
    maxnorm = fmax(maxnorm, creal(norm[i]));

  double mine = 0.0;
  for (i=0; i<d; i++)
    mine = fmin(mine, creal(e[i]));

  // merits

  struct merit {
    double val;
    int idx;
  } merits[d];

  // how should we sort the eigenstates ?
  // according to energy, according to overlap with intrinsic state ?
  // problematic as ordering may change during minimization

  for (i=0; i<d; i++) {
    merits[i].idx = i;
    if (creal(norm[i]) < Min.minnormkmix*maxnorm)
      merits[i].val = -100000.0;
    else
      merits[i].val = e[i]/mine + 0.1*norm[i]/maxnorm;
  }

  // sort according to merits value
  qsort(merits, d, sizeof(struct merit), cmpmerit);

  // now write good eigenvalues and eigenvectors into egood, Vgood

  i=0; while(i<d && merits[i].val > -100000.0) {
    egood[i] = e[merits[i].idx];
    for (k=0; k<j+1; k++)
      Vgood[k+i*(j+1)] = V[k+merits[i].idx*(j+1)];
    i++;
  }
  *ngood = i;

}

/*
void kmixingstates(complex double* H, complex double* N, int j, double thresh, 
                   double* egood, complex double* Vgood, int* ngood)
{
  complex double e[j+1];
  complex double V[(j+1)*(j+1)];
  int d;

  generalizedeigensystem(H, N, j+1, thresh, e, V, &d);

  // we got d eigenstates, filter out the noise
  // see sortEigenstates

  double norm[d];
  int i, k, l;

  for (i=0; i<d; i++) {
      norm[i] = 0.0;
      for (l=0; l<j+1; l++)
        for (k=0; k<j+1; k++)
          norm[i] += conj(V[k+i*(j+1)])*N[k+l*(j+1)]*V[l+i*(j+1)];
    }

  double maxnorm = 0.0;
  for (i=0; i<d; i++)
    maxnorm = fmax(maxnorm, creal(norm[i]));

  // merits

  struct merit {
    double val;
    int idx;
  } merits[d];

  for (i=0; i<d; i++) {
    merits[i].idx = i;
    if (creal(norm[i]) < Min.minnormkmix*maxnorm)
      merits[i].val = -100000.0;
    else
      merits[i].val = -creal(e[i]);
  }

  // sort according to merits value
  qsort(merits, d, sizeof(struct merit), cmpmerit);

  // now write good eigenvalues and eigenvectors into egood, Vgood

  i=0; while(i<d && merits[i].val > -100000.0) {
    egood[i] = -merits[i].val;
    for (k=0; k<j+1; k++)
      Vgood[k+i*(j+1)] = V[k+merits[i].idx*(j+1)];
    i++;
  }
  *ngood = i;

}
*/

#ifndef MPI
void calcprojectedHamiltonianmatrix(const SlaterDet* Q, const SlaterDet* Qp, 
                                    const Interaction* Int,
                                    int j, int par, 
                                    angintegrationpara* angpara,
                                    complex double* H, complex double* N, complex double* J2)
{
  SlaterDet* Qpp = Work.Qpp;
  SlaterDetAux* X = Work.X;

  int nang = angpara->n;

  reinitangintegration(Q, Qp, 0, 0, angpara);

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      J2[idxjmk(j,k,kp)] = 0.0;
    }

  complex double nme, hme, j2me;
  double alpha, beta, gamma, w;

  int i, ip;
  for (i=0; i<nang; i++) {
    copySlaterDet(Qp, Qpp);

    getangintegrationpoint(i, angpara, &alpha, &beta, &gamma, &w);

    rotateSlaterDet(Qpp, alpha, beta, gamma);
    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);
      
      calcSlaterDetAuxod(Q, Qpp, X);
      nme = X->ovlap;
      calcHamiltonianod(Int, Q, Qpp, X, &hme);
      calcConstraintJ2od(Q, Qpp, X, &j2me);

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

}
#endif

#ifdef MPI
void calcprojectedHamiltonianmatrix(const SlaterDet* Q, const SlaterDet* Qp,
                                    const Interaction* Int,
                                    int j, int par, 
                                    angintegrationpara* angpara,
                                    complex double* H, complex double* N, complex double* J2)
{
  int nang = angpara->n;

  reinitangintegration(Q, Qp, 0, 0, angpara);

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      J2[idxjmk(j,k,kp)] = 0.0;
    }
   
  int task = TASKHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Q);
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

}
#endif


#ifndef MPI
void calcgradprojectedHamiltonianmatrix(const SlaterDet* Q, const SlaterDet* Qp, 
                                        const Interaction* Int,
                                        int j, int par, 
                                        angintegrationpara* angpara,
                                        complex double* H, complex double* N,
                                        gradSlaterDet* dH, gradSlaterDet* dN)
{
  SlaterDet* Qpp = Work.Qpp;
  SlaterDetAux* X = Work.X;
  gradSlaterDet* dh = Work.dhproj;
  gradSlaterDet* dn = Work.dnproj;
  gradSlaterDetAux* dX = Work.dX;

  int nang = angpara->n;

  reinitangintegration(Q, Qp, 0, 0, angpara);

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      zerogradSlaterDet(&dH[idxjmk(j,k,kp)]);
      zerogradSlaterDet(&dN[idxjmk(j,k,kp)]);
    }

  double alpha, beta, gamma, w;

  int i, ip;
  for (i=0; i<nang; i++) {
    copySlaterDet(Qp, Qpp);

    getangintegrationpoint(i, angpara, &alpha, &beta, &gamma, &w);

    rotateSlaterDet(Qpp, alpha, beta, gamma);

    for (ip=0; ip<=1; ip++) {
      if (ip) invertSlaterDet(Qpp);

      calcSlaterDetAuxod(Q, Qpp, X);
      calcgradSlaterDetAuxod(Q, Qpp, X, dX);
      calcgradOvlapod(Q, Qpp, X, dX, dn);
      dn->val = X->ovlap;
      calcgradHamiltonianod(Int, Q, Qpp, X, dX, dh);

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

}
#endif

#ifdef MPI
void calcgradprojectedHamiltonianmatrix(const SlaterDet* Q, const SlaterDet* Qp, 
                                        const Interaction* Int,
                                        int j, int par, 
                                        angintegrationpara* angpara,
                                        complex double* H, complex double* N,
                                        gradSlaterDet* dH, gradSlaterDet* dN)
{
  gradSlaterDet* dh = Work.dhproj;
  gradSlaterDet* dn = Work.dnproj;

  int nang = angpara->n;

  reinitangintegration(Q, Qp, 0, 0, angpara);

  int k,kp;
  for (kp=-j; kp<=j; kp=kp+2)
    for (k=-j; k<=j; k=k+2) {
      H[idxjmk(j,k,kp)] = 0.0; 
      N[idxjmk(j,k,kp)] = 0.0; 
      zerogradSlaterDet(&dH[idxjmk(j,k,kp)]);
      zerogradSlaterDet(&dN[idxjmk(j,k,kp)]);
    }

  int task = TASKGRADHAMILTONIANOD;
  BroadcastTask(&task);

  BroadcastSlaterDet(Q);
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
}
#endif


void initWork(int j, int par,
              angintegrationpara* AngPar,
              const Interaction* Int,
              const SlaterDet* Q)
{       
  int k;

  Work.Qp = malloc(sizeof(SlaterDet));
  Work.Qpp = malloc(sizeof(SlaterDet));
  Work.X = malloc(sizeof(SlaterDetAux));
  Work.dX = malloc(sizeof(gradSlaterDetAux));

  initSlaterDet(Q, Work.Qp);
  initSlaterDet(Q, Work.Qpp);
  initSlaterDetAux(Q, Work.X);
  initgradSlaterDetAux(Q, Work.dX);

  Work.dH = malloc((j+1)*(j+1)*sizeof(gradSlaterDet));
  Work.dN = malloc((j+1)*(j+1)*sizeof(gradSlaterDet));
    
  for (k=0; k<(j+1)*(j+1); k++) {
    initgradSlaterDet(Q, &Work.dH[k]); 
    initgradSlaterDet(Q, &Work.dN[k]);
  }

  Work.dhproj = malloc(sizeof(gradSlaterDet));
  Work.dnproj = malloc(sizeof(gradSlaterDet));

  initgradSlaterDet(Q, Work.dhproj);
  initgradSlaterDet(Q, Work.dnproj);

  Work.V = malloc((j+1)*(j+1)*sizeof(complex double));

  Work.H = malloc((j+1)*(j+1)*sizeof(complex double));
  Work.N = malloc((j+1)*(j+1)*sizeof(complex double));
  Work.J2 = malloc((j+1)*(j+1)*sizeof(complex double));

}

// perform multiconfiguration calculation with fixed intrinsic states
// return eigenvalue #ival

// assumes fix matrix elements to be available in Work

void calcprojectedHamiltonian(const SlaterDet* Q,
                              const Interaction* Int,
                              int j, int par,
                              int ival,
                              angintegrationpara* angpara,
                              double* eintr, double* eproj)
{
  complex double* H = Work.H;
  complex double* N = Work.N;
  complex double* J2 = Work.J2;
  complex double* V = Work.V;
  int dim = Work.dim;

  // Q is intrinsic state to be varied
  // orient properly
  // then calculate needed matrix elements
  // solve eigenvalue problem in smaller subspace spanned by K-mixing eigenstates

  SlaterDet* Qp = Work.Qp;
  SlaterDetAux* X = Work.X;

#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif
   
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

  // intrinsic state

  calcSlaterDetAuxod(Qp, Qp, X);

  if (isnan(X->ovlap)) {
    *eintr = ENERGYPENALTY;
    *eproj = ENERGYPENALTY;
    return;
  }

  nintr = X->ovlap;
  calcHamiltonianod(Int, Qp, Qp, X, &hintr);

  // projected matrix elements for Qp
  calcprojectedHamiltonianmatrix(Qp, Qp, Int, j, par, angpara,
                                 H, N, J2);

  // K-mixing for Qp
  double ed[j+1];

  kmixingstates(H, N, j, Min.threshkmix, ed, V, &dim);

  if (!dim) {
    fprintf(stderr, "no good eigenstate for Q !\n");
  }

  complex double N2[(j+1)*(j+1)];
  multcmat(N, N, N2, j+1);

  complex double hproj[j+1], nproj[j+1], n2proj[j+1], j2proj[j+1];

  int i;
  for (i=0; i<dim; i++) {
    hproj[i] = expect(j+1, H, &V[i*(j+1)]);
    nproj[i] = expect(j+1, N, &V[i*(j+1)]);
    n2proj[i] = expect(j+1, N2, &V[i*(j+1)]);
    j2proj[i] = expect(j+1, J2, &V[i*(j+1)]);
  }

  for (i=0; i<dim; i++)
    fprintf(stderr, "Q eigenstate %d: E = %8.3f MeV, n = %8.7f, J2 = %6.3f\n", i, 
            hbc*creal(hproj[i]/nproj[i]), creal(n2proj[i]/(nintr*nproj[i])), creal(j2proj[i]/nproj[i]));

  *eintr = hintr/nintr;
  *eproj = creal(hproj[ival]/nproj[ival]);
}



// gradient will be written to Work.dhproj

void calcgradprojectedHamiltonian(const SlaterDet* Q, 
				  const Interaction* Int,
				  int j, int par, 
                                  int ival,
				  angintegrationpara* angpara,
				  double* eintr, double* eproj)
{
  SlaterDet* Qp = Work.Qp;
  SlaterDetAux* X = Work.X;
  complex double* H = Work.H;
  complex double* N = Work.N;
  complex double* J2 = Work.J2;
  complex double* V = Work.V;
  int dim = Work.dim;
  gradSlaterDet* dH = Work.dH;
  gradSlaterDet* dN = Work.dN;
  gradSlaterDet* dhproj = Work.dhproj;
  gradSlaterDet* dnproj = Work.dnproj;

#ifdef ORIENTED
  double alpha0, beta0, gamma0;
#endif

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

  // projected matrix elements for Qp
  calcgradprojectedHamiltonianmatrix(Qp, Qp, Int, j, par, angpara,
                                     H, N,
                                     dH, dN);

  // K-mixing for Qp
  double ed[j+1];

  kmixingstates(H, N, j, Min.threshkmix,
                ed, V, &dim);

  if (!dim) {
    fprintf(stderr, "no good eigenstate for Q !\n");
  }

  complex double N2[(j+1)*(j+1)];
  multcmat(N, N, N2, j+1);

  complex double hproj[j+1], nproj[j+1], n2proj[j+1], j2proj[j+1];

  int i;
  for (i=0; i<dim; i++) {
    hproj[i] = expect(j+1, H, &V[i*(j+1)]);
    nproj[i] = expect(j+1, N, &V[i*(j+1)]);
    n2proj[i] = expect(j+1, N2, &V[i*(j+1)]);
    j2proj[i] = expect(j+1, J2, &V[i*(j+1)]);
  }

  for (i=0; i<dim; i++)
    fprintf(stderr, "Q eigenstate %d: E = %8.3f MeV, n = %8.7f, J2 = %6.3f\n", i, 
            hbc*creal(hproj[i]/nproj[i]), creal(n2proj[i]/(nintr*nproj[i])), creal(j2proj[i]/nproj[i]));

  *eproj = creal(hproj[ival]/nproj[ival]);

  // gradients

  gradexpect(j+1, dH, &V[ival*(j+1)], dhproj);
  gradexpect(j+1, dN, &V[ival*(j+1)], dnproj);

  multgradSlaterDet(dhproj, 1.0/nproj[ival]);
  addmulttogradSlaterDet(dhproj, dnproj, -hproj[ival]/(nproj[ival]*nproj[ival]));

#ifdef ORIENTED
  rotategradSlaterDet(dhproj, alpha0, beta0, gamma0);
#endif
}				  



static void copyxtopara(const double* x, Para* q)
{
  int i;
  for (i=0; i<q->n; i++)
    q->x[i] = x[i];
}


void MinimizeDONLP2vapp(const Interaction* Int, 
                        int j, int par, int ival,
                        double threshkmix, double minnormkmix,
                        double alpha,
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
  Min.ival = ival;
  Min.threshkmix = threshkmix;
  Min.minnormkmix = minnormkmix;
  Min.alpha = alpha;
  Min.projpar = projpar;
  Min.Const = Const;
  Min.nconst = nconst;
  Min.q = q;
  Min.P = P;

  Min.Q = malloc(sizeof(SlaterDet));
  P->ParainitSlaterDet(q, Min.Q); 

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

  // nail to center-of-mass
  moveboostSlaterDet(Min.Q, Work.X);

  double eintr, eproj;

  calcprojectedHamiltonian(Min.Q,
			   Min.Int,
			   Min.j, Min.par, Min.ival,
                           Min.projpar,
                           &eintr, &eproj);

  *fx = (eproj+Min.alpha*eintr);
  FORTRAN(o8cnt).icf++;

  fprintf(stderr, "\tme      E = % 8.3f MeV, Eproj = %8.3f MeV, Eintr = %8.3f MeV,\n", 
	  hbc*(*fx), hbc*eproj, hbc*eintr);

}


void FORTRAN(egradf)(const double* x, double* gradf)
{
#ifndef MPI
  if (sigterminate)
    longjmp(env, 1);
#endif

  copyxtopara(x, Min.q);
  Min.P->ParatoSlaterDet(Min.q, Min.Q);

  // nail to center-of-mass
  moveboostSlaterDet(Min.Q, Work.X);

  double eintr, eproj;

  calcgradprojectedHamiltonian(Min.Q,
			       Min.Int,
			       Min.j, Min.par, Min.ival, 
                               Min.projpar,
			       &eintr, &eproj);

  // add gradient from intrinsic energy

  calcSlaterDetAux(Min.Q, Work.X);
  calcgradSlaterDetAux(Min.Q, Work.X, Work.dX);
  calcgradHamiltonian(Min.Int, Min.Q, Work.X, Work.dX, Work.dH);

  addmulttogradSlaterDet(Work.dhproj, Work.dH, Min.alpha);

  Min.P->ParaprojectgradSlaterDet(Min.q, Work.dhproj, gradf);
  FORTRAN(o8cnt).icgf++;
  fprintf(stderr, "grad %3d: \tE = %8.3f MeV, Eproj = %8.3f MeV, Eintr = %8.3f MeV\n", 
	  FORTRAN(o8cnt).icgf, hbc*(eproj+Min.alpha*eintr), hbc*eproj, hbc*eintr);


  if (Min.log && !(FORTRAN(o8cnt).icgf % Min.log)) {
    char logfname[255];
    sprintf(logfname, "%s.minvapp.%03d.log", Min.logfile, FORTRAN(o8cnt).icgf);
    FILE* logfp;
    if (!(logfp = fopen(logfname, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", logfname);
    } else {
      fprintinfo(logfp);
      fprintf(logfp, "# step %3d: \tE = %8.3f MeV, Eproj = %8.3f MeV, Eintr = %8.3f MeV\n", 
	      FORTRAN(o8cnt).icgf, hbc*(eproj+Min.alpha*eintr), hbc*eproj, hbc*eintr);
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

  if (*i > 1) {
    // nail to center-of-mass
    moveboostSlaterDet(Min.Q, Work.X);
  }

  calcSlaterDetAux(Min.Q, Work.X);
  Min.Const[*i-1].me(Min.Q, Work.X, hxi);

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

  if (*i > 1) {
    // nail to center-of-mass
    moveboostSlaterDet(Min.Q, Work.X);
  }

  calcSlaterDetAux(Min.Q, Work.X);
  calcgradSlaterDetAux(Min.Q, Work.X, Work.dX);
  zerogradSlaterDet(Work.dH);
  Min.Const[*i-1].gradme(Min.Q, Work.X, Work.dX, Work.dH);

  Min.P->ParaprojectgradSlaterDet(Min.q, Work.dH, gradhi);
  fprintf(stderr, "\t\t\tgrad %4s = %8.3f\n", 
	  Min.Const[*i-1].label, Min.Const[*i-1].output(Work.dH->val));
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
