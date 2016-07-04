/**

   \file OneNucleonOvlaps.c

   calculate spectroscopic amplitudes in momentum space


   (c) 2004-2012 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"

#include "OneNucleonOvlaps.h"

#include "numerics/cmath.h"
#include "numerics/cmat.h"
#include "numerics/gaussquad.h"
#include "numerics/sphericalharmonics.h"
#include "numerics/clebsch.h"
#include "numerics/lapack.h"

#include "misc/utils.h"
#include "misc/physics.h"

#define SQR(x)	((x)*(x))


char* speclabel[NSPEC] = { "S12", "P12", "P32", "D32", "D52" };
int   specj[NSPEC]     = { 1, 1, 3, 3, 5 };
int   specl[NSPEC]     = { 0, 2, 2, 4, 4 };
int   io[NSPEC]        = { 0, 2, 4, 8, 12 };
int   dimo[NSPEC]      = { 2, 2, 4, 4, 6 };


// The following Many-Body Operators have to initialized with
// int dim and FormfactorPara par before use

// Array of Multipole Formfactor Operators
ManyBodyOperator OpOneNucleonOvlap[NSPEC];


ManyBodyOperators OpOneNucleonOvlaps = {
  n : NSPEC,
  Op : OpOneNucleonOvlap,
  me : calcOneNucleonOvlapsod
};



void initOpOneNucleonOvlaps(OneNucleonOvlapsPara* par)
{
  char* name = malloc(60);
  sprintf(name, "OneNucleonOvlaps-%05.2f-%d-%d-%d%s",
	  par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");

  OpOneNucleonOvlaps.name = name;
  OpOneNucleonOvlaps.dim = par->npoints;
  OpOneNucleonOvlaps.size = par->npoints;
  OpOneNucleonOvlaps.par = par;

  int i;
  for (i=0; i<NSPEC; i++) {
    char* name = malloc(60);
    sprintf(name, "OneNucleonOvlap%s-%05.2f-%d-%d-%d%s",
	    speclabel[i],
	    par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");
    OpOneNucleonOvlap[i].name = name;
    OpOneNucleonOvlap[i].rank = specj[i];
    OpOneNucleonOvlap[i].pi = specl[i] % 4 ? 1 : 0;
    OpOneNucleonOvlap[i].dim = par->npoints;
    OpOneNucleonOvlap[i].size = par->npoints;
    OpOneNucleonOvlap[i].par = par;
    OpOneNucleonOvlap[i].me = NULL;
  }
}


// static SlaterDet and overlap matrices
static SlaterDet *QAp;
static complex double *n, *nlu;


// overlap between Gaussians
static complex double calcGaussianOvlap(const Gaussian* G1, const Gaussian* G2)
{
  if (G1->xi != G2->xi)
    return 0.0;

  int T;
  complex double lambda, alpha, pi[3], pi2, S, R;
    
  lambda = 1.0/(conj(G1->a)+G2->a);
  alpha = conj(G1->a)*G2->a*lambda;
  for (int i=0; i<3; i++)	
    pi[i] = lambda*I*(conj(G1->b[i]) - G2->b[i]);
  pi2 = cvec3sqr(pi);
  
  T = 1;
  S = conj(G1->chi[0])*G2->chi[0] + conj(G1->chi[1])*G2->chi[1];
  R = cpow32(2*M_PI*alpha)*cexp(0.5*pi2/lambda);
  return (T*S*R);
}


// calculate overlap matrix between QA and QB
static void calcOvlapsAB(const SlaterDet* QA, const SlaterDet* QB, complex double* N)
{
  assert(QA->A+1 == QB->A);
  
  int A=QB->A;
  int* idxA=QA->idx; int* idxB=QB->idx; 
  int* ngA=QA->ng; int* ngB=QB->ng;
  Gaussian* GA=QA->G; Gaussian* GB=QB->G;

  int k,l,ki,li;

  // overlap matrix between QB and QA

  for (l=0; l<A; l++)
    for (k=0; k<A-1; k++) {
      N[k+l*A] = 0.0;
      if (GA[idxA[k]].xi == GB[idxB[l]].xi) { 
	for (li=0; li<ngB[l]; li++)
	  for (ki=0; ki<ngA[k]; ki++) {
	    N[k+l*A] += calcGaussianOvlap(&GA[idxA[k]+ki], &GB[idxB[l]+li]); 
	  }
      }
    }
}


// calculate overlap matrix for QA+nucleon with QB

static void calcSpecOvlaps(const SlaterDet* QA,
		    int iso, double q[3],
		    const SlaterDet* QB,
		    complex double* N, complex double* Nlu,
		    complex double spec[2])
{
  assert(QA->A+1 == QB->A);
  
  int A=QB->A;
  int ngaussB=QB->ngauss; 
  int* idxB=QB->idx; int* ngB=QB->ng;
  Gaussian* GB=QB->G;

  int l,li;

  double q2 = vec3sqr(q);

  int ipiv[A];
  int info;

  // overlap matrix elements of QB with QA should have been calculated already

  // overlap matrix elements with plane wave
  // isospin overlap zero or one by definition

  // spatial overlap is independent from spin: calculate only once
  complex double nR[ngaussB];
  Gaussian* g;

  for (li=0; li<ngaussB; li++) {
    g = &GB[li];
    if (iso != g->xi)
      nR[li] = 0.0;
    else
      nR[li] = cpow32(g->a)* cexp(-0.5*g->a*q2 - I*(q[0]*g->b[0]+q[1]*g->b[1]+q[2]*g->b[2]));
  }

  // be careful: chi[0] is the spin up, chi[1] the spin-down component
  // but destroying a spin-up is a spin-down operator
  // tensor operator matrix elements go from 0: -m ... m

  // spin-up plane wave in bra
  for (l=0; l<A; l++) {
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[0]* nR[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[1]);
  
  // spin-down plane wave in bra
  for (l=0; l<A; l++) {
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[1]* nR[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[0]);

}


// hermitian adjoint of destruction operator is tricky
// tilde{a}_j,m = (-1)^(j+m) a_j,-m

void calcOneNucleonOvlapsod(OneNucleonOvlapsPara* par,
			    const SlaterDet* QA, const SlaterDet* QB,
			    const SlaterDetAux* dummyX,
			    complex double* specamplitude)
{
  double qmax = par->qmax;
  int ialpha, nalpha = par->nalpha;
  int icosb, ncosb = par->ncosb;
  int npoints = par->npoints;
  int m, ml, ms, o, i;
  int j, l;
  double alpha, beta, weight;
  double q, k[3];
  complex double sa[2];
  int gaussoverlapsdone = 0;

  // initialize workspace if necessary
  if (!QAp) {
    QAp = malloc(sizeof(SlaterDet));
    initSlaterDet(QA, QAp);
    
    n = malloc(SQR(QB->A)*sizeof(complex double));
    nlu = malloc(SQR(QB->A)*sizeof(complex double));
  }

  copySlaterDet(QA, QAp);

  // add single nucleon with momentum q
  // proton or neutron ?
  int iso = (QA->Z == QB->Z) ? -1 : +1;
  double massA = QA->Z*mproton+QA->N*mneutron;

  for (o=0; o<NSPEC; o++)
    for (i=0; i<dimo[o]*npoints; i++)
      specamplitude[i+io[o]*npoints] = 0.0;

  double alphal[nalpha], walphal[nalpha];
  double cosbl[ncosb], wcosbl[ncosb];
  
  ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0.0, alphal, walphal);
  GaussLegendrePoints(ncosb, -1, 1, cosbl, wcosbl); 

  for (i=0; i<npoints; i++) {
    q = qmax* i/(npoints-1);
    for (ialpha=0; ialpha<nalpha; ialpha++)
      for (icosb=0; icosb<ncosb; icosb++) {
        alpha = alphal[ialpha];
        beta = acos(cosbl[icosb]);
        weight = walphal[ialpha]*wcosbl[icosb];
 
        k[0] = q*cos(alpha)*sin(beta); 
        k[1] = q*sin(alpha)*sin(beta);
        k[2] = q*cos(beta);

	if (par->recoil) {
	  // boost QA by -k
	  double VcmA[3] = {-k[0]/massA, -k[1]/massA, -k[2]/massA};
	  copySlaterDet(QA, QAp);
	  boostSlaterDet(QAp, VcmA);
	  gaussoverlapsdone = 0;
	}

	if (!gaussoverlapsdone) {
	  calcOvlapsAB(QAp, QB, n);
	  gaussoverlapsdone = 1;
	}

	calcSpecOvlaps(QAp, iso, k, QB, n, nlu, sa);

        for (ms=-1; ms<=1; ms=ms+2) {
          for (o=0; o<NSPEC; o++) {
            j = specj[o]; l = specl[o];
            for (m=-j; m<=j; m=m+2) {
              ml = m-ms;
              if (-l<=ml && ml<=l) {
                specamplitude[(-m+j)/2 + i*dimo[o] + io[o]*npoints] += weight* 
                  ((j-m)%4 ? -1 : +1)* 
                  clebsch(l, 1, j, ml, ms, m)* sa[(ms+1)/2]* 
                  conj(Y(l,ml,beta,alpha));
              }
            }
          }
        }
    }	
  }
}



// output routine


void writeOneNucleonOvlaps(FILE* fp,
			   const Projection* P,
			   const OneNucleonOvlapsPara* p,
			   int jfin, int pfin, int afin,
			   int jini, int pini, int aini,
			   void* specamplitudeme,
			   const Eigenstates* Efin,
			   const Eigenstates* Eini)
{
  int jmax = P->jmax;
  int ipjfin = idxpij(jmax,pfin,jfin);
  int ipjini = idxpij(jmax,pini,jini);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = Efin->index[ipjfin][afin];
  int idxini = Eini->index[ipjini][aini];
  double normfin = Efin->norm[ipjfin][idxfin];
  double normini = Eini->norm[ipjini][idxini];
  complex double (****specampme)[npoints] = specamplitudeme;
  complex double *sa = specampme[ipjfin][ipjini][idxfin][idxini];
  int i;
  double q;
  complex double ml;
  
  // calculate spectroscopic factor
  double dq = qmax/(npoints-1);
  double S = 0.0;
  for (i=0; i<npoints; i++) {
    q = i*dq;
    S += dq*SQR(q)* (jfin+1.0)/(jini+1.0)*SQR(cabs(sa[i]))/(normfin*normini);
  }
  fprintf(stderr, "Spectroscopic factor: %6.3f\n", S);

  // write spectroscopic amplitudes
  for (i=0; i<npoints; i++) {
    q = i*dq;
    
    // we include factor 1/sqrt(2*Jini+1) into definiton
    // of spectroscopic amplitude

    ml = sqrt((jfin+1.0)/(jini+1.0))* sa[i]/sqrt(normfin*normini);

    fprintf(fp, "%6.3f\t%13.5g\t%13.5g\n", q, creal(ml), cimag(ml));
  }
}
