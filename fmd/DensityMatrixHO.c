/**

   \file DensityMatrixHO.c

   calculate density matrix in harmonic oscillator basis


   (c) 2006 Thomas Neff

*/


#include <stdlib.h>
#include <complex.h>

#include "SlaterDet.h"
#include "HOBasis.h"
#include "DensityMatrixHO.h"


// workspace

complex double ampg1[4*HODIMMAX];
complex double ampg2[4*HODIMMAX];


// we are recacalculating the amplitudes for the single-particle
// states again and again !

// not ease to solve in present framework, let's hope that it's
// fast enough like it is

static void ob_densmatho(DensityMatrixHOPar* par,
			 const Gaussian* G1, const Gaussian* G2,
			 const GaussianAux* X, complex double* rho)
{
  int nmax = par->nmax;
  double omega = par->omega;
  int dim = par->dim;

  amplitudesHOGaussian(G1, nmax, omega, ampg1);
  amplitudesHOGaussian(G2, nmax, omega, ampg2);
  
  int i,j;
  for (j=0; j<dim; j++)
   for (i=0; i<dim; i++)
     rho[i+j*dim] += conj(ampg1[i])*ampg2[j];
}
  

// density matrix needs not be real ?
void calcDensityMatrixHO(const DensityMatrixHOPar* par, 
			 const SlaterDet* Q, const SlaterDetAux* X,
			 complex double* rho)
{
  int d = par->dim;

  OneBodyOperator op_ob_densmatho = { dim: d*d, opt: 0, par: par, 
				      me: ob_densmatho };

  calcSlaterDetAuxod(Q, Q, X);
  calcSlaterDetOBMEod(Q, Q, X, &op_ob_densmatho, rho);

  int i;
  for (i=0; i<d*d; i++)
    rho[i] /= X->ovlap;
}


/*
void calcDensityMatrixHOod(const DensityMatrixHOPar* par,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   const SlaterDetAux* X,
			   complex double* rho)
{
  int d = par->dim;

  OneBodyOperator op_ob_densmatho = { dim: d*d, opt: 0, par: par, 
				      me: ob_densmatho };

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_densmatho, rho);
}
*/


 
// workspace

int ngwork = 0;
complex double *ampgwork;
complex double *ampgpwork;


void calcDensityMatrixHOod(const DensityMatrixHOPar* par,
			   const SlaterDet* Q, const SlaterDet* Qp, 
			   const SlaterDetAux* X,
			   complex double* rho)
{
  int A=Q->A; 
  int ngauss=Q->ngauss; int ngaussp=Qp->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  complex double ovl = X->ovlap;

  int nmax = par->nmax;
  double omega = par->omega;
  int dim = par->dim;

  if (ngwork != Q->A*MAXNG) {
    ngwork = Q->A*MAXNG;
    ampgwork = malloc(ngwork*dim*sizeof(complex double));
    ampgpwork = malloc(ngwork*dim*sizeof(complex double));
  }

  complex double (*ampg)[dim] = ampgwork;
  complex double (*ampgp)[dim] = ampgpwork;

  int i;
  for (i=0; i<ngauss; i++)
    amplitudesHOGaussian(&G[i], nmax, omega, ampg[i]);
  for (i=0; i<ngaussp; i++)
    amplitudesHOGaussian(&Gp[i], nmax, omega, ampgp[i]);

  int k,l, ki,li, ik,il;

  for (i=0; i<dim*dim; i++)
    rho[i] = 0.0;

  for (l=0; l<A; l++)
    for (k=0; k<A; k++) {

      for (li=0; li<ngp[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  
	  for (il=0; il<dim; il++)
	    for (ik=0; ik<dim; ik++)
	      rho[ik+il*dim] += 
		conj(ampg[idx[k]+ki][ik])* ampgp[idx[l]+li][il]*
		o[l+k*A]*ovl;

    }		
}	
