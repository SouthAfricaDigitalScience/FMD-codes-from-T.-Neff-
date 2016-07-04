/**

   \file TwoNucleonOvlaps.c

   calculate two-nucleon spectroscopic amplitudes in momentum space


   (c) 2012 Thomas Neff

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

#include "TwoNucleonOvlaps.h"

#include "numerics/cmath.h"
#include "numerics/cmat.h"
#include "numerics/gaussquad.h"
#include "numerics/sphericalharmonics.h"
#include "numerics/clebsch.h"
#include "numerics/lapack.h"

#include "misc/utils.h"
#include "misc/physics.h"

#define SQR(x)	((x)*(x))


// only T=1 amplitudes implemented

// JLlSj

char* specTlabel[NSPECT] = { "0Ss00", "0Pp11", "0Dd02", "0Ff13", "0Gg04" };
int   specTJ[NSPECT]     = { 0, 0, 0, 0, 0 };
int   specTL[NSPECT]     = { 0, 2, 4, 6, 8 };
int   specTl[NSPECT]     = { 0, 2, 4, 6, 8 };
int   specTS[NSPECT]     = { 0, 2, 0, 2, 0 };
int   specTj[NSPECT]     = { 0, 2, 4, 6, 8 };
int   ioT[NSPECT]        = { 0, 1, 2, 3, 4 };
int   dimoT[NSPECT]      = { 1, 1, 1, 1, 1 };

// Jj1j2

char* specYlabel[NSPECY] = { "0s12s12", "0p12p12", "0p32p32", "0d32d32", "0d52d52" };
int   specYJ[NSPECY]     = { 0, 0, 0, 0, 0 };
int   specYj1[NSPECY]    = { 1, 1, 3, 3, 5 };
int   specYl1[NSPECY]    = { 0, 2, 2, 4, 4 };
int   specYj2[NSPECY]    = { 1, 1, 3, 3, 5 };
int   specYl2[NSPECY]    = { 0, 2, 2, 4, 4 };
int   ioY[NSPECY]        = { 0, 1, 2, 3, 4 };
int   dimoY[NSPECY]      = { 1, 1, 1, 1, 1 };


// The following Many-Body Operators have to initialized with
// int dim and FormfactorPara par before use

// Array of Multipole Formfactor Operators
ManyBodyOperator OpTwoNucleonOvlapT[NSPECT];

ManyBodyOperator OpTwoNucleonOvlapY[NSPECY];


ManyBodyOperators OpTwoNucleonOvlapsT = {
  n : NSPECT,
  Op : OpTwoNucleonOvlapT,
  me : calcTwoNucleonOvlapsTod
};

ManyBodyOperators OpTwoNucleonOvlapsY = {
  n : NSPECY,
  Op : OpTwoNucleonOvlapY,
  me : calcTwoNucleonOvlapsYod
};


void initOpTwoNucleonOvlapsT(TwoNucleonOvlapsPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoNucleonOvlapsT-%05.2f-%d-%d-%d%s",
	  par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");

  OpTwoNucleonOvlapsT.name = name;
  OpTwoNucleonOvlapsT.dim = SQR(par->npoints);
  OpTwoNucleonOvlapsT.size = SQR(par->npoints);
  OpTwoNucleonOvlapsT.par = par;

  int i;
  for (i=0; i<NSPECT; i++) {
    char* name = malloc(60);
    sprintf(name, "TwoNucleonOvlapsT%s-%05.2f-%d-%d-%d%s",
	    specTlabel[i],
	    par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");
    OpTwoNucleonOvlapT[i].name = name;
    OpTwoNucleonOvlapT[i].rank = specTJ[i];
    OpTwoNucleonOvlapT[i].pi = (specTL[i]+specTl[i]) % 4 ? 1 : 0;
    OpTwoNucleonOvlapT[i].dim = SQR(par->npoints);
    OpTwoNucleonOvlapT[i].size = SQR(par->npoints);
    OpTwoNucleonOvlapT[i].par = par;
    OpTwoNucleonOvlapT[i].me = NULL;
  }
}

void initOpTwoNucleonOvlapsY(TwoNucleonOvlapsPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoNucleonOvlapsY-%05.2f-%d-%d-%d%s",
	  par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");

  OpTwoNucleonOvlapsY.name = name;
  OpTwoNucleonOvlapsY.dim = SQR(par->npoints);
  OpTwoNucleonOvlapsY.size = SQR(par->npoints);
  OpTwoNucleonOvlapsY.par = par;

  int i;
  for (i=0; i<NSPECY; i++) {
    char* name = malloc(60);
    sprintf(name, "TwoNucleonOvlapsY%s-%05.2f-%d-%d-%d%s",
	    specYlabel[i],
	    par->qmax, par->npoints, par->nalpha, par->ncosb, par->recoil ? "-recoil" : "");
    OpTwoNucleonOvlapY[i].name = name;
    OpTwoNucleonOvlapY[i].rank = specYJ[i];
    OpTwoNucleonOvlapY[i].pi = (specYl1[i]+specYl2[i]) % 4 ? 1 : 0;
    OpTwoNucleonOvlapY[i].dim = SQR(par->npoints);
    OpTwoNucleonOvlapY[i].size = SQR(par->npoints);
    OpTwoNucleonOvlapY[i].par = par;
    OpTwoNucleonOvlapY[i].me = NULL;
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
  assert(QA->A+2 == QB->A);
  
  int A=QB->A;
  int* idxA=QA->idx; int* idxB=QB->idx; 
  int* ngA=QA->ng; int* ngB=QB->ng;
  Gaussian* GA=QA->G; Gaussian* GB=QB->G;

  int k,l,ki,li;

  // overlap matrix between QB and QA

  for (l=0; l<A; l++)
    for (k=0; k<A-2; k++) {
      N[k+l*A] = 0.0;
      if (GA[idxA[k]].xi == GB[idxB[l]].xi) { 
	for (li=0; li<ngB[l]; li++)
	  for (ki=0; ki<ngA[k]; ki++) {
	    N[k+l*A] += calcGaussianOvlap(&GA[idxA[k]+ki], &GB[idxB[l]+li]); 
	  }
      }
    }
}


// calculate overlap matrix for QA+2nucleons with QB

static void calcSpecOvlaps(const SlaterDet* QA,
			   int iso1, double q1[3],
			   int iso2, double q2[3],
			   const SlaterDet* QB,
			   complex double* N, complex double* Nlu,
			   complex double spec[2][2])
{
  assert(QA->A+2 == QB->A);
  
  int A=QB->A;
  int ngaussB=QB->ngauss; 
  int* idxB=QB->idx; int* ngB=QB->ng;
  Gaussian* GB=QB->G;

  int l,li;

  double q12 = vec3sqr(q1);
  double q22 = vec3sqr(q2);

  int ipiv[A];
  int info;

  // overlap matrix elements of QB with QA should have been calculated already

  // overlap matrix elements with plane wave
  // isospin overlap zero or one by definition

  // spatial overlap is independent from spin: calculate only once
  complex double nR1[ngaussB], nR2[ngaussB];
  Gaussian* g;

  for (li=0; li<ngaussB; li++) {
    g = &GB[li];
    if (iso1 != g->xi)
      nR1[li] = 0.0;
    else
      nR1[li] = cpow32(g->a)* cexp(-0.5*g->a*q12 - I*(q1[0]*g->b[0]+q1[1]*g->b[1]+q1[2]*g->b[2]));
    if (iso2 != g->xi)
      nR2[li] = 0.0;
    else
      nR2[li] = cpow32(g->a)* cexp(-0.5*g->a*q22 - I*(q2[0]*g->b[0]+q2[1]*g->b[1]+q2[2]*g->b[2]));
  }

  // be careful: chi[0] is the spin up, chi[1] the spin-down component
  // but destroying a spin-up is a spin-down operator
  // tensor operator matrix elements go from 0: -m ... m

  // spin-up/spin-up plane waves in bra
  for (l=0; l<A; l++) {
    N[(A-2)+l*A] = 0.0;
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-2)+l*A] += GB[idxB[l]+li].chi[0]* nR1[idxB[l]+li];
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[0]* nR2[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[1][1]);

  // spin-up/spin-down plane waves in bra
  for (l=0; l<A; l++) {
    N[(A-2)+l*A] = 0.0;
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-2)+l*A] += GB[idxB[l]+li].chi[0]* nR1[idxB[l]+li];
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[1]* nR2[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[1][0]);

  // spin-down/spin-up plane waves in bra
  for (l=0; l<A; l++) {
    N[(A-2)+l*A] = 0.0;
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-2)+l*A] += GB[idxB[l]+li].chi[1]* nR1[idxB[l]+li];
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[0]* nR2[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[0][1]);
  
  // spin-down/spin-down plane waves in bra
  for (l=0; l<A; l++) {
    N[(A-2)+l*A] = 0.0;
    N[(A-1)+l*A] = 0.0;
    for (li=0; li<ngB[l]; li++) {
      N[(A-2)+l*A] += GB[idxB[l]+li].chi[1]* nR1[idxB[l]+li];
      N[(A-1)+l*A] += GB[idxB[l]+li].chi[1]* nR2[idxB[l]+li];
    }
  }

  copycmat(A, N, Nlu);
  FORTRAN(zgetrf)(&A, &A, Nlu, &A, ipiv, &info);
  FORTRAN(zdet)(Nlu, &A, &A, ipiv, &spec[0][0]);

}


// hermitian adjoint of destruction operator is tricky
// tilde{a}_j,m = (-1)^(j+m) a_j,-m

// to be done
void calcTwoNucleonOvlapsTod(TwoNucleonOvlapsPara* par,
			     const SlaterDet* QA, const SlaterDet* QB,
			     const SlaterDetAux* dummyX,
			     complex double* specamplitude)
{
  double qmax = par->qmax;
  int ialphaQ, ialphaq, nalpha = par->nalpha;
  int icosbQ, icosbq, ncosb = par->ncosb;
  int npoints = par->npoints; int npoints2 = SQR(npoints);
  int ML, ml, ms1, ms2, MS, mj, M, o, iQ, iq;
  int L, l, S, j, J;
  double alphaQ, alphaq, betaQ, betaq, weightQ, weightq;
  double Ql, ql, Q[3], q[3], k1[3], k2[3];
  complex double sa[2][2];
  int gaussoverlapsdone = 0;

  // initialize workspace if necessary
  if (!QAp) {
    QAp = malloc(sizeof(SlaterDet));
    initSlaterDet(QA, QAp);
    
    n = malloc(SQR(QB->A)*sizeof(complex double));
    nlu = malloc(SQR(QB->A)*sizeof(complex double));
  }

  copySlaterDet(QA, QAp);

  // add two nucleons with momenta q1 and q2
  // proton/proton or neutron/neutron
  int iso1, iso2;
  if (QA->Z == QB->Z) {
    iso1 = -1; iso2 = -1;
  } else if (QA->Z+1 == QB->Z) {
    fprintf(stderr, "pn case not implemented yet !\n");
    exit(-1);
  } else if (QA->Z+2 == QB->Z) {
    iso1 = +1; iso2 = +1;
  }

  double massA = QA->Z*mproton+QA->N*mneutron;
  double mass1 = (iso1 == +1) ? mproton : mneutron;
  double mass2 = (iso2 == +1) ? mproton : mneutron;

  for (o=0; o<NSPECT; o++)
    for (int i=0; i<dimoT[o]*npoints2; i++)
      specamplitude[i+ioT[o]*npoints2] = 0.0;

  double alphal[nalpha], walphal[nalpha];
  double cosbl[ncosb], wcosbl[ncosb];
  
  ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0.0, alphal, walphal);
  GaussLegendrePoints(ncosb, -1, 1, cosbl, wcosbl); 

  for (iQ=0; iQ<npoints; iQ++) {
    Ql = 2*qmax* iQ/(npoints-1);
    for (ialphaQ=0; ialphaQ<nalpha; ialphaQ++)
      for (icosbQ=0; icosbQ<ncosb; icosbQ++) {
	alphaQ = alphal[ialphaQ];
	betaQ = acos(cosbl[icosbQ]);
	weightQ = walphal[ialphaQ]*wcosbl[icosbQ];

	Q[0] = Ql*cos(alphaQ)*sin(betaQ); 
	Q[1] = Ql*sin(alphaQ)*sin(betaQ);
	Q[2] = Ql*cos(betaQ);

	if (par->recoil) {
	  // boost QA by -Q
	  double VcmA[3] = {-Q[0]/massA, -Q[1]/massA, -Q[2]/massA};
	  copySlaterDet(QA, QAp);
	  boostSlaterDet(QAp, VcmA);
	  gaussoverlapsdone = 0;
	}

	if (!gaussoverlapsdone) {
	  calcOvlapsAB(QAp, QB, n);
	  gaussoverlapsdone = 1;
	}

	for (iq=0; iq<npoints; iq++) {
	  ql = qmax* iq/(npoints-1);
	  for (ialphaq=0; ialphaq<nalpha; ialphaq++)
	    for (icosbq=0; icosbq<ncosb; icosbq++) {
	      alphaq = alphal[ialphaq];
	      betaq = acos(cosbl[icosbq]);
	      weightq = walphal[ialphaq]*wcosbl[icosbq];

	      q[0] = ql*cos(alphaq)*sin(betaq); 
	      q[1] = ql*sin(alphaq)*sin(betaq);
	      q[2] = ql*cos(betaq);

	      for (int i=0; i<3; i++) {
		k1[i] = q[i] + mass1/(mass1+mass2)*Q[i];
		k2[i] = -q[i] + mass2/(mass1+mass2)*Q[i];
	      }

	      calcSpecOvlaps(QAp, iso1, k1, iso2, k2, QB, n, nlu, sa);

	      for (ms1=-1; ms1<=1; ms1=ms1+2) {
		for (ms2=-1; ms2<=1; ms2=ms2+2) {

		  for (o=0; o<NSPECT; o++) {
		    l = specTl[o]; S = specTS[o]; j = specTj[o];
		    L = specTL[o]; J = specTJ[o];
		    
		    MS = ms1+ms2;
		    if (-S<=MS && MS <= S) {
		      for (M=-J; M<=J; M=M+2) {
			for (ML=-L; ML<=L; ML=ML+2) {
			  mj = M-ML;
			  if (-j<=mj && mj<=j) {
			    ml = mj-MS;
			    if (-l<=ml && ml<=l) {

			      // are the phase factors correct ?

			      specamplitude[(-M+J)/2 + (iq+iQ*npoints)*dimoY[o] + ioY[o]*npoints2] += 1.0/sqrt(2.0)*weightQ*weightq*
				((1-ms1)%4 ? -1 : +1)* ((1-ms2)%4 ? -1 : +1)*
				((l-ml)%4 ? -1 : +1)* ((L-ML)%4 ? -1 : +1)*
				clebsch(L, j, J, ML, mj, M)*
				clebsch(l, S, j, ml, MS, mj)*
				clebsch(1, 1, S, ms1, ms2, MS)* 
				sa[(ms1+1)/2][(ms2+1)/2]* 
				conj(Y(L,ML,betaQ,alphaQ))*
				conj(Y(l,ml,betaq,alphaq));
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
            }
	}
      }
  }	
}


void calcTwoNucleonOvlapsYod(TwoNucleonOvlapsPara* par,
			     const SlaterDet* QA, const SlaterDet* QB,
			     const SlaterDetAux* dummyX,
			     complex double* specamplitude)
{
  double qmax = par->qmax;
  int ialpha1, ialpha2, nalpha = par->nalpha;
  int icosb1, icosb2, ncosb = par->ncosb;
  int npoints = par->npoints; int npoints2 = SQR(npoints);
  int m1, m2, ml1, ml2, ms1, ms2, M, o, i1, i2;
  int j1, j2, l1, l2, J;
  double alpha1, alpha2, beta1, beta2, weight1, weight2;
  double q1l, q2l, q1[3], q2[3], k1[3], k2[3];
  complex double sa[2][2];
  int gaussoverlapsdone = 0;

  // initialize workspace if necessary
  if (!QAp) {
    QAp = malloc(sizeof(SlaterDet));
    initSlaterDet(QA, QAp);
    
    n = malloc(SQR(QB->A)*sizeof(complex double));
    nlu = malloc(SQR(QB->A)*sizeof(complex double));
  }

  copySlaterDet(QA, QAp);

  // add two nucleons with momenta q1 and q2
  // proton/proton or neutron/neutron
  int iso1, iso2;
  if (QA->Z == QB->Z) {
    iso1 = -1; iso2 = -1;
  } else if (QA->Z+1 == QB->Z) {
    fprintf(stderr, "pn case not implemented yet !\n");
    exit(-1);
  } else if (QA->Z+2 == QB->Z) {
    iso1 = +1; iso2 = +1;
  }

  double massA = QA->Z*mproton+QA->N*mneutron;
  double massB = QB->Z*mproton+QB->N*mneutron;
  double mass1 = (iso1 == +1) ? mproton : mneutron;
  double mass2 = (iso2 == +1) ? mproton : mneutron;

  for (o=0; o<NSPECY; o++)
    for (int i=0; i<dimoY[o]*npoints2; i++)
      specamplitude[i+ioY[o]*npoints2] = 0.0;

  double alphal[nalpha], walphal[nalpha];
  double cosbl[ncosb], wcosbl[ncosb];
  
  ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0.0, alphal, walphal);
  GaussLegendrePoints(ncosb, -1, 1, cosbl, wcosbl); 

  for (i1=0; i1<npoints; i1++) {
    q1l = qmax* i1/(npoints-1);
    for (i2=0; i2<npoints; i2++) {
      q2l = qmax* i2/(npoints-1);

      for (ialpha1=0; ialpha1<nalpha; ialpha1++)
	for (icosb1=0; icosb1<ncosb; icosb1++) {
	  alpha1 = alphal[ialpha1];
	  beta1 = acos(cosbl[icosb1]);
	  weight1 = walphal[ialpha1]*wcosbl[icosb1];
	  
	  q1[0] = q1l*cos(alpha1)*sin(beta1); 
	  q1[1] = q1l*sin(alpha1)*sin(beta1);
	  q1[2] = q1l*cos(beta1); 

	  for (ialpha2=0; ialpha2<nalpha; ialpha2++)
	    for (icosb2=0; icosb2<ncosb; icosb2++) {
	      alpha2 = alphal[ialpha2];
	      beta2 = acos(cosbl[icosb2]);
	      weight2 = walphal[ialpha2]*wcosbl[icosb2];

	      q2[0] = q2l*cos(alpha2)*sin(beta2); 
	      q2[1] = q2l*sin(alpha2)*sin(beta2);
	      q2[2] = q2l*cos(beta2);

	      if (par->recoil) {
		for (int i=0; i<3; i++) {
		  k1[i] = q1[i] - mass1/(mass1+massA)*q2[i];
		  k2[i] = q2[i];
		} 
	      } else {
		for (int i=0; i<3; i++) {
		  k1[i] = q1[i];
		  k2[i] = q2[i];
		}
	      }

	      if (par->recoil) {
		// boost QA by -(k1+k2)
		double VcmA[3] = {-(k1[0]+k2[0])/massA, -(k1[1]+k2[1])/massA, -(k1[2]+k2[2])/massA};
		copySlaterDet(QA, QAp);
		boostSlaterDet(QAp, VcmA);
		gaussoverlapsdone = 0;
	      }

	      if (!gaussoverlapsdone) {
		calcOvlapsAB(QAp, QB, n);
		gaussoverlapsdone = 1;
	      }

	      calcSpecOvlaps(QAp, iso1, k1, iso2, k2, QB, n, nlu, sa);

	      for (ms1=-1; ms1<=1; ms1=ms1+2) {
		for (ms2=-1; ms2<=1; ms2=ms2+2) {
		  for (o=0; o<NSPECY; o++) {
		    j1 = specYj1[o]; l1 = specYl1[o];
		    j2 = specYj2[o]; l2 = specYl2[o];
		    J = specYJ[o];
		    
		    for (M=-J; M<=J; M=M+2) {
		      for (m1=-j1; m1<=j1; m1=m1+2) {
			m2 = M-m1;
			if (-j2<=m2 && m2<=j2) {
			  ml1 = m1-ms1;
			  if (-l1<=ml1 && ml1<=l1) {
			    ml2 = m2-ms2;
			    if (-l2<=ml2 && ml2<=l2) {

			      specamplitude[(-M+J)/2 + (i1+i2*npoints)*dimoY[o] + ioY[o]*npoints2] += 1.0/sqrt(2.0)*weight1*weight2*
				clebsch(j1, j2, J, m1, m2, M)*
				((j1-m1)%4 ? -1 : +1)* ((j2-m2)%4 ? -1 : +1)*
				clebsch(l1, 1, j1, ml1, ms1, m1)*
				clebsch(l2, 1, j2, ml2, ms2, m2)* 
				sa[(ms1+1)/2][(ms2+1)/2]* 
				conj(Y(l1,ml1,beta1,alpha1))*
				conj(Y(l2,ml2,beta2,alpha2));
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
            }
	}
    }
  }	
}

// output routine


void writeTwoNucleonOvlapsT(FILE* fp,
			    const Projection* P,
			    const TwoNucleonOvlapsPara* p,
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
  int iQ, iq;
  double Q, q;
  complex double ml;
  
  // calculate spectroscopic factor
  double dQ = 2*qmax/(npoints-1);
  double dq = qmax/(npoints-1);

  double S = 0.0;
  for (iQ=0; iQ<npoints; iQ++) {
    Q = iQ*dQ;
    for (iq=0; iq<npoints; iq++) {
      q = iq*dq;

      S += dq*SQR(q)*dQ*SQR(Q)* (jfin+1.0)/(jini+1.0)*
	SQR(cabs(sa[iq+iQ*npoints]))/(normfin*normini);
    }
  }
  fprintf(stderr, "Spectroscopic factor: %6.3f\n", S);

  for (iQ=0; iQ<npoints; iQ++) {
    Q = iQ*dQ;

    for (iq=0; iq<npoints; iq++) {
      q = iq*dq;

      ml = sqrt((jfin+1.0)/(jini+1.0))* sa[iq+iQ*npoints]/sqrt(normfin*normini);

      fprintf(fp, "(%13.5g, %13.5g)\t", creal(ml), cimag(ml));
    }
    fprintf(fp, "\n");
  }

}


void writeTwoNucleonOvlapsY(FILE* fp,
			    const Projection* P,
			    const TwoNucleonOvlapsPara* p,
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
  int i1, i2;
  double q1, q2;
  complex double ml;
  
  // calculate spectroscopic factor
  double dq = qmax/(npoints-1);

  double S = 0.0;
  for (i2=0; i2<npoints; i2++) {
    q2 = i2*dq;
    for (i1=0; i1<npoints; i1++) {
      q1 = i1*dq;

      S += dq*SQR(q1)*dq*SQR(q2)* (jfin+1.0)/(jini+1.0)*
	SQR(cabs(sa[i1+i2*npoints]))/(normfin*normini);
    }
  }
  fprintf(stderr, "Spectroscopic factor: %6.3f\n", S);

  for (i1=0; i1<npoints; i1++) {
    q1 = i1*dq;

    for (i2=0; i2<npoints; i2++) {
      q2 = i2*dq;

      ml = sqrt((jfin+1.0)/(jini+1.0))* sa[i1+i2*npoints]/sqrt(normfin*normini);

      fprintf(fp, "(%13.5g, %13.5g)\t", creal(ml), cimag(ml));
    }
    fprintf(fp, "\n");
  }

}
