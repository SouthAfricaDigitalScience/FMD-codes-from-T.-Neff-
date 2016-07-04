/**

   \file ProjectedShellOccupations.c

   calculate the many-body occupations in a harmonic oscillator basis
   see Suzuki et al, PRC 54, 2073 (1996)

   in principle we would like to calculate < | exp {I theta H^HO_intr} | >,
   as H^HO_intr is a two-body operator that is not possible, approximate treatment
   by dividing with <Qcm|exp{I theta H^HO}|Qcm'> assuming factorization of |Q> into
   |Qintr> and |Qcm>
   should be exact for AMD type wave functions

   proton/neutron occupation numbers not correct yet !
   We could calculate occupations with respect to proton/neutron center-of-mass

   (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"
#include "Projection.h"
#include "ProjectedShellOccupations.h"

#include "misc/physics.h"
#include "numerics/cmath.h"
#include "numerics/cmat.h"
#include "numerics/lapack.h"


ManyBodyOperator OpShellOccupations = {
 name : NULL,
 rank : 0,
 pi : 0,
 dim : 0,
 size : 0,
 par : NULL,
 me : calcShellOccupationsod
};

void initOpShellOccupations(ShellOccupationsPara* par)
{
  char* shelloccupname = malloc(40);
  sprintf(shelloccupname, "ShellOccupations-%05.2f-%d-%d",
	  hbc*par->omega, par->nmax, par->nint);
  OpShellOccupations.name = shelloccupname;
  OpShellOccupations.par = par;
  OpShellOccupations.dim = 3*(par->nmax+1);
  OpShellOccupations.size = 3*(par->nmax+1);
}


// determinant of matrix A, matrix gets destroyed
static complex double cdet(int n, complex double* A)
{
  int N = n;
  int ipiv[N];
  int info;
  complex double det;
  
  FORTRAN(zgetrf)(&N, &N, A, &N, ipiv, &info);
  FORTRAN(zdet)(A, &N, &N, ipiv, &det);

  return det;
}


complex double calcovl(complex double a1, const complex double* b1,
                       complex double a2, const complex double* b2)
{
  complex double lambda, alpha; 
  complex double pi[3];

  lambda = 1.0/(conj(a1)+a2);
  alpha = lambda* conj(a1)*a2;
  for (int i=0; i<3; i++)
    pi[i] = I*lambda*(conj(b1[i])-b2[i]);

  return cpow32(2*M_PI*alpha)*cexp(0.5/lambda*cvec3sqr(pi));
}


complex double calcb(complex double a1, const complex double* b1,
                     complex double a2, const complex double* b2,
                     double aho, complex double z)
{
  // might be more elegant to give the full matrix element
  // <a1,b1| exp{i theta (HHO-3/2)} |a2,b2>

  // now do it by:
  // aplying exp(HHO) on right hand side |a2,b2>
  // calculate overlap with left hand side <a1,b1|

  complex double a2p = aho*(aho+a2-(aho-a2)*z*z)/(aho+a2+(aho-a2)*z*z);
  complex double b2p[3];
  for (int i=0; i<3; i++)
    b2p[i] = 2*aho*z/(aho+a2+(aho-a2)*z*z)* b2[i];
  
  complex double lambdap = 1/(conj(a1)+a2p);
  complex double alphap = conj(a1)*a2p*lambdap;
  complex double pip[3];
  for (int i=0; i<3; i++)
    pip[i] = I*lambdap*(conj(b1[i])-b2p[i]);
  
  return cpow32(2*a2/(aho+a2-(aho-a2)*z*z))*
    cexp(-0.5*(1-z*z)/(aho+a2+(aho-a2)*z*z)*cvec3sqr(b2))*
    cpow32(2*M_PI*alphap)*cexp(0.5*cvec3sqr(pip)/lambdap);
}


void calcbqq(const Gaussian* G1, const Gaussian* G2,
	     const GaussianAux* X,
	     double omega, complex double z,
	     complex double *b, complex double *bp, complex double *bn)
{
  if (!X->T) {
    *b = 0.0; *bp = 0.0; *bn = 0.0;
    return;
  }

  double aho = 1.0/(mass(G1->xi)*omega);
  complex double me;
  
  me = calcb(G1->a, G1->b, G2->a, G2->b, aho, z);

  *b = me*X->S*X->T;
    
  if (G1->xi == -1)
    *bp = X->R* X->S* X->T;
  else
    *bp = me* X->S* X->T;

  if (G1->xi == +1)
    *bn = X->R* X->S* X->T;
  else
    *bn = me* X->S* X->T;

}

// 0hbw oscillator quanta

int nq[] = {  0,                                 // no nucleons
              0,0,                               // s-shell
              1,2,3,4,5,6,                       // p-shell
              8,10,12,14,16,18,20,22,24,26,      // sd-shell
              29,32,35,38,41,44,47,50,53,56,     // pf-shell
              59,62,65,68,71,74,77,80,83,86 };


static complex double* B;
static complex double* Bp;
static complex double* Bn;


void calcShellOccupationsod(ShellOccupationsPara* par,
			    const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    complex double *shelloccupations)
{
  int A = Q->A; int ngauss = Q->ngauss;
  int *idx = Q->idx, *idxp = Qp->idx;
  int *ng = Q->ng, *ngp = Qp->ng;
  Gaussian *G = Q->G, *Gp = Qp->G;
  GaussianAux *Gaux = X->Gaux;

  int nmax = par->nmax;
  int nint = par->nint;
  double omega = par->omega;

  int n0p = nq[Q->Z];
  int n0n = nq[Q->N];
  int n0 = n0p+n0n;

  // correct for excitations of CM wave functions
  // procedure is exact if intrinsic and CM wave functions factorize

  double M, Mp, Mn;
  double tcm, tcmp;
  double acm, acmp;
  double xcm[3], xcmp[3];
  double vcm[3], vcmp[3];
  complex double bcm[3], bcmp[3];

  M = mproton*Q->Z+mneutron*Q->N;
  Mp = mproton*Q->Z;
  Mn = mneutron*Q->N;

  SlaterDetAux Xx;
  initSlaterDetAux(Q, &Xx);

  calcSlaterDetAux(Q, &Xx);
  calcTCM(Q, &Xx, &tcm);
  calcCMPosition(Q, &Xx, xcm);
  calcCMVelocity(Q, &Xx, vcm);

  calcSlaterDetAux(Qp, &Xx);
  calcTCM(Qp, &Xx, &tcmp);
  calcCMPosition(Qp, &Xx, xcmp);
  calcCMVelocity(Qp, &Xx, vcmp);

  freeSlaterDetAux(&Xx);

  int i;
  acm = 0.75/(tcm*M);
  for (i=0; i<3; i++)
    bcm[i] = xcm[i]+I*acm*M*vcm[i];

  acmp = 0.75/(tcmp*M);
  for (i=0; i<3; i++)
    bcmp[i] = xcmp[i]+I*acmp*M*vcmp[i];

  double ahocm, ahocmp, ahocmn; 
  ahocm = 1.0/(M*omega);
  ahocmp = 1.0/(Mp*omega);
  ahocmn = 1.0/(Mn*omega);

  complex double ovlcm;
  ovlcm = calcovl(acm, bcm, acmp, bcmp);

  complex double Bcm, Bcmp, Bcmn;
  
  complex double (*occupations)[nmax+1] = shelloccupations;

  // B matrix workspace

  if (!B) {
    B = initcmat(A);
    Bp = initcmat(A);
    Bn = initcmat(A);
  }

  int n, k, l, ki, li;
  int it;
  double theta;
  complex double z;
  complex double b, bp, bn;

  for (it=0; it<3; it++) 
    for (n=0; n<=nmax; n++)
      occupations[it][n] = 0.0;

  for (i=0; i<nint; i++) {
    theta = i*2*M_PI/nint;
    z = cexp(I*theta);

    // approximate contribution from cm HO Hamiltonian
    Bcm = calcb(acm, bcm, acmp, bcmp, ahocm, z);
    Bcmp = calcb(acm, bcm, acmp, bcmp, ahocmp, z);
    Bcmn = calcb(acm, bcm, acmp, bcmp, ahocmn, z);

    // calculate B matrices

    for (l=0; l<A; l++)
      for (k=0; k<A; k++) {
    
	B[k+l*A] = 0.0; Bp[k+l*A] = 0.0; Bn[k+l*A] = 0.0;

	for (li=0; li<ngp[l]; li++)
	  for (ki=0; ki<ng[k]; ki++) {

	    calcbqq(&G[idx[k]+ki], &Gp[idxp[l]+li], 
		    &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss],
		    omega, z,
		    &b, &bp, &bn);

	    B[k+l*A] += b; Bp[k+l*A] += bp; Bn[k+l*A] += bn;
	  }
      }

    // determinants of B matrices

    complex double detB = cdet(A, B);
    complex double detBp = cdet(A, Bp);
    complex double detBn = cdet(A, Bn);

    for (n=0; n<=nmax; n++) {
      occupations[0][n] += 1.0/nint* cexp(-I*theta*(n+n0))* ovlcm/Bcm* detB;
      occupations[1][n] += 1.0/nint* cexp(-I*theta*(n+n0p))* ovlcm/Bcm* detBp;
      occupations[2][n] += 1.0/nint* cexp(-I*theta*(n+n0n))* ovlcm/Bcm* detBn;
    }

  }
}

static double dabs(double x)
{
  return (x<0 ? -x : x);
}

static void fprintpercent(FILE* fp, double x)
{
  // if (dabs(100*x) < 0.01)
  //   fprintf(fp, "          *");
  // else
  fprintf(fp, "   %8.2f", 100*x);
}  

void writeprojectedShellOccupations(FILE* fp, const ShellOccupationsPara* par,
				    const Projection* P,
				    void *shelloccupations,
				    const Eigenstates* E,
				    int totalonly)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int nmax=par->nmax;

  int p,j,i,n;

  char prefix[8];

  int* idx;
  int ngood;
  complex double *norm, *H;
  complex double (**occupations)[3*(nmax+1)] = shelloccupations;
  complex double (*occ)[3*(nmax+1)];
  

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

	if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
	else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
	norm = E->norm[idxpij(jmax,p,j)];
	H = E->v[idxpij(jmax,p,j)];
	occ = occupations[idxpij(jmax,p,j)];

	fprintf(fp, "\n%s N      = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
	fprintf(fp, "\n%s H      = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));

	double sum[ngood], mean[ngood], occu;

	for (i=0; i<ngood; i++) {
	  sum[i] = 0.0; mean[i] = 0.0;
	}
	for (n=0; n<=nmax; n++) {
	  fprintf(fp, "\n%s Q[%02d]  = ", prefix, n); 
	  for (i=0; i<ngood; i++) {
	    occu = creal(occ[idx[i]][n]/norm[idx[i]]);
	    sum[i] += occu;
	    mean[i] += n*occu;
	    fprintpercent(fp, occu);
	  }
	}
	fprintf(fp, "\n%s Q[--]  = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintpercent(fp, sum[i]);
	fprintf(fp, "\n%s Qmean  = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", mean[i]);
	
	if (!totalonly) {

	  for (i=0; i<ngood; i++) {
	    sum[i] = 0.0;
	    mean[i] = 0.0;
	  }
	  for (n=0; n<=nmax; n++) {
	    fprintf(fp, "\n%s Qp[%02d] = ", prefix, n); 
	    for (i=0; i<ngood; i++) {
	      occu = creal(occ[idx[i]][n+(nmax+1)]/norm[idx[i]]);
	      sum[i] += occu;
	      mean[i] += n*occu;
	      fprintpercent(fp, occu);
	    }
	  }
	  fprintf(fp, "\n%s Qp[--] = ", prefix); 
	  for (i=0; i<ngood; i++) 
	    fprintpercent(fp, sum[i]);
	  fprintf(fp, "\n%s Qpmean = ", prefix); 
	  for (i=0; i<ngood; i++) 
	    fprintf(fp, "   %8.3f", mean[i]);

	  for (i=0; i<ngood; i++) {
	    sum[i] = 0.0;
	    mean[i] = 0.0;
	  }
	  for (n=0; n<=nmax; n++) {
	    fprintf(fp, "\n%s Qn[%02d] = ", prefix, n); 
	    for (i=0; i<ngood; i++) {
	      occu = creal(occ[idx[i]][n+2*(nmax+1)]/norm[idx[i]]);
	      sum[i] += occu;
	      mean[i] += n*occu;
	      fprintpercent(fp, occu);
	    }
	  }
	  fprintf(fp, "\n%s Qn[--] = ", prefix); 
	  for (i=0; i<ngood; i++) 
	    fprintpercent(fp, sum[i]);
	  fprintf(fp, "\n%s Qnmean = ", prefix); 
	  for (i=0; i<ngood; i++) 
	    fprintf(fp, "   %8.3f", mean[i]);
	}

	fprintf(fp, "\n");

      }
    }
}
