/**

  \file gradSlaterDet.c

  gradients for SlaterDet matrixelements


  (c) 2003 Thomas Neff

*/

#include <assert.h>
#include <stdlib.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "gradGaussian.h"
#include "gradSlaterDet.h"

#include "numerics/cmath.h"

#define SQR(x) (x)*(x)


void allocategradSlaterDet(gradSlaterDet* G, int A)
{
  G->gradval = malloc(MAXNG*A*sizeof(gradGaussian));
}

void initgradSlaterDet(const SlaterDet* Q, gradSlaterDet* G)
{
  G->ngauss = Q->ngauss;
  G->gradval = (gradGaussian*) malloc(MAXNG*Q->A*sizeof(gradGaussian));
  zerogradSlaterDet(G);
}

void freegradSlaterDet(gradSlaterDet* G)
{
  free(G->gradval);
}


void zerogradSlaterDet(gradSlaterDet* G)
{
  int i;

  G->val = 0.0;
  for (i=0; i<G->ngauss; i++)
    zerogradGaussian(&G->gradval[i]);
}


void copygradSlaterDet(const gradSlaterDet* G, gradSlaterDet* Gp)
{
  assert(G->ngauss == Gp->ngauss);

  int i;

  Gp->val = G->val;
  for (i=0; i<G->ngauss; i++)
    Gp->gradval[i] = G->gradval[i];
}


void multgradSlaterDet(gradSlaterDet* G, complex double x)
{
  int i;

  G->val *= x;
  for (i=0; i<G->ngauss; i++)
    multgradGaussian(&G->gradval[i], x);
}


void addtogradSlaterDet(gradSlaterDet* G, const gradSlaterDet* dG)
{
  int i;

  G->val += dG->val;
  for (i=0; i<G->ngauss; i++)
    addtogradGaussian(&G->gradval[i], &dG->gradval[i]);
}


void addmulttogradSlaterDet(gradSlaterDet* G, const gradSlaterDet* dG, 
			    complex double x)
{
  int i;

  G->val += x*dG->val;
  for (i=0; i<G->ngauss; i++)
    addmulttogradGaussian(&G->gradval[i], &dG->gradval[i], x);
}


void rotategradSlaterDet(gradSlaterDet* G, double alpha, double beta, double gamma)
{
  double euler[3] = { alpha, beta, gamma };

  complex double R2[2][2];
  rotatemat2(euler, R2);

  double R3[3][3];
  rotatemat3(euler, R3);
  
  int i;
  for (i=0; i<G->ngauss; i++)
    rotategradGaussian(&G->gradval[i], R3, R2);

}

void allocategradSlaterDetAux(gradSlaterDetAux* dX, int A)
{
  dX->dGaux = malloc(SQR(MAXNG*A)*sizeof(gradGaussianAux));
  dX->dno = malloc(MAXNG*SQR(A)*sizeof(gradGaussian));
}


void initgradSlaterDetAux(const SlaterDet* Q, gradSlaterDetAux* dX)
{
  dX->ngauss = Q->ngauss;
  dX->dGaux = malloc(SQR(MAXNG*Q->A)*sizeof(gradGaussianAux));
  dX->dno = malloc(MAXNG*SQR(Q->A)*sizeof(gradGaussian)); 
}  
  

void freegradSlaterDetAux(gradSlaterDetAux* dX)
{
  free(dX->dGaux);
  free(dX->dno);
}


void calcgradSlaterDetAux(const SlaterDet* Q, const SlaterDetAux* X,
			  gradSlaterDetAux* dX)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;

  int k,l,m,mi;
  gradGaussian dn;

  for (l=0; l<ngauss; l++)
    for (k=0; k<ngauss; k++)
      calcgradGaussianAux(&G[k], &G[l], &Gaux[k+l*ngauss],
			  &dGaux[k+l*ngauss]);

  for (l=0; l<A; l++)
    for (k=0; k<ngauss; k++)
      zerogradGaussian(&dno[k+l*ngauss]);

  for (m=0; m<A; m++)
    for (k=0; k<ngauss; k++) {
      zerogradGaussian(&dn);
      for (mi=0; mi<ng[m]; mi++)
	addtogradGaussian(&dn, &dGaux[k+(idx[m]+mi)*ngauss].dQ);
      for (l=0; l<A; l++)
	addmulttogradGaussian(&dno[k+l*ngauss], &dn, o[m+l*A]);
    }
}


void calcgradSlaterDetAuxod(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X,
			    gradSlaterDetAux* dX)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* idxp=Qp->idx;
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;

  int k,ki,l,li,m,mi;
  gradGaussian dn;

  for (l=0; l<A; l++)
    for (li=0; li<ngp[l]; li++)
      for (k=0; k<A; k++)
	for (ki=0; ki<ng[k]; ki++)
	  calcgradGaussianAux(&G[idx[k]+ki], &Gp[idxp[l]+li], 
			      &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss],
			      &dGaux[(idx[k]+ki)+(idxp[l]+li)*ngauss]);

  for (l=0; l<A; l++)
    for (k=0; k<ngauss; k++)
      zerogradGaussian(&dno[k+l*ngauss]);

  for (m=0; m<A; m++)
    for (k=0; k<ngauss; k++) {
      zerogradGaussian(&dn);
      for (mi=0; mi<ngp[m]; mi++)
	addtogradGaussian(&dn, &dGaux[k+(idxp[m]+mi)*ngauss].dQ);
      for (l=0; l<A; l++)
	addmulttogradGaussian(&dno[k+l*ngauss], &dn, o[m+l*A]);
    }
}


void calcgradOvlapod(const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, const gradSlaterDetAux* dX,
		     gradSlaterDet* grad)
{
  int A=Q->A; int ngauss=Q->ngauss;
  int* idx = Q->idx; int* ng = Q->ng;
  gradGaussian* dno = dX->dno;
  gradGaussian* dval = grad->gradval;
  complex double ovl = X->ovlap;

  int k,ki;

  grad->ngauss = Q->ngauss;
  zerogradSlaterDet(grad);

  for (k=0; k<A; k++)
    for (ki=0; ki<ng[k]; ki++)
      addmulttogradGaussian(&dval[idx[k]+ki], &dno[(idx[k]+ki)+k*ngauss], 
			    ovl);

  grad->val = ovl;
}


void calcgradSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   const gradOneBodyOperator* op, 
			   gradSlaterDet* grad)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* ng = Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double* val = &grad->val; gradGaussian* dval = grad->gradval;

  int k,l,s,ki,li,si;
  complex double gval;
  gradGaussian gdval;

  for (l=0; l<A; l++)
    for (k=0; k<A; k++)
      if (!op->opt || Gaux[idx[k]+idx[l]*ngauss].T) {
	gval = 0.0;
	for (ki=0; ki<ng[k]; ki++) {
	  zerogradGaussian(&gdval);
	  for (li=0; li<ng[l]; li++) {
	    op->me(op->par, 
		   &G[idx[k]+ki], &G[idx[l]+li], 
		   &Gaux[(idx[k]+ki)+(idx[l]+li)*ngauss],
		   &dGaux[(idx[k]+ki)+(idx[l]+li)*ngauss],
		   &gval, &gdval);
	  }
	  addmulttogradGaussian(&dval[idx[k]+ki], &gdval, o[l+k*A]); 
	}
	*val += gval*o[l+k*A];
	for (s=0; s<A; s++)
	  for (si=0; si<ng[s]; si++)
	    addmulttogradGaussian(&dval[idx[s]+si], &dno[(idx[s]+si)+k*ngauss],
				  -gval*o[l+s*A]);
      }
}


void calcgradSlaterDetOBMEod(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const gradOneBodyOperator* op, 
			     gradSlaterDet* grad)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* idxp = Qp->idx;
  int* ng = Q->ng; int* ngp = Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double ovl = X->ovlap;
  complex double val; gradGaussian* dval = grad->gradval;

  int k,l,s,ki,li,si;
  complex double gval;
  gradGaussian gdval;

  val = 0.0;
  for (l=0; l<A; l++)
    for (k=0; k<A; k++)
      if (!op->opt || Gaux[idx[k]+idxp[l]*ngauss].T) {
	gval = 0.0;
	for (ki=0; ki<ng[k]; ki++) {
	  zerogradGaussian(&gdval);
	  for (li=0; li<ngp[l]; li++) {
	    op->me(op->par, 
		   &G[idx[k]+ki], &Gp[idxp[l]+li], 
		   &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss],
		   &dGaux[(idx[k]+ki)+(idxp[l]+li)*ngauss],
		   &gval, &gdval);
	  }
	  addmulttogradGaussian(&dval[idx[k]+ki], &gdval, o[l+k*A]*ovl); 
	}
	val += gval*o[l+k*A]*ovl;
	for (s=0; s<A; s++)
	  for (si=0; si<ng[s]; si++)
	    addmulttogradGaussian(&dval[idx[s]+si], &dno[(idx[s]+si)+k*ngauss],
				  -gval*o[l+s*A]*ovl);
      }

  grad->val += val;

  // this term is comming from the derivative of the overlap
  for (k=0; k<A; k++)
    for (ki=0; ki<ng[k]; ki++)
      addmulttogradGaussian(&dval[idx[k]+ki], &dno[(idx[k]+ki)+k*ngauss], val);
}


void calcgradSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   const gradTwoBodyOperator* op, 
			   gradSlaterDet* grad)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* ng = Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double* val = &grad->val; gradGaussian* dval = grad->gradval;

  int k,l,m,n,s, ki,li,mi,ni,si;
  complex double gval;
  gradGaussian gdval;
  complex double ooa;

  for (n=0; n<A; n++)
    for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idx[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idx[m]*ngauss].T) {

	      ooa = o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A];

	      gval = 0.0;

	      for (ki=0; ki<ng[k]; ki++) {
		zerogradGaussian(&gdval);
		for (ni=0; ni<ng[n]; ni++)
		  for (li=0; li<ng[l]; li++)
		    for (mi=0; mi<ng[m]; mi++)

		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &G[idx[m]+mi], &G[idx[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idx[n]+ni)*ngauss],
			     &dGaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &gval, &gdval);

		addmulttogradGaussian(&dval[idx[k]+ki], &gdval, ooa);
	      }
	      *val += 0.5*gval*ooa;

	      for (s=0; s<A; s++) {
		ooa = o[m+s*A]*o[n+l*A]-o[n+s*A]*o[m+l*A];
		for (si=0; si<ng[s]; si++)
		  addmulttogradGaussian(&dval[idx[s]+si], 
					&dno[(idx[s]+si)+k*ngauss],
					-gval*ooa);
	      }
	    }
      }
}


void calcgradSlaterDetTBMEod(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const gradTwoBodyOperator* op, 
			     gradSlaterDet* grad)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* idxp = Qp->idx; 
  int* ng = Q->ng; int* ngp = Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double val; gradGaussian* dval = grad->gradval;
  complex double ovl = X->ovlap;

  int k,l,m,n,s, ki,li,mi,ni,si;
  complex double gval;
  gradGaussian gdval;
  complex double ooa;

  val = 0.0;
  for (n=0; n<A; n++)
    for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idxp[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idxp[m]*ngauss].T) {

	      ooa = o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A];

	      gval = 0.0;

	      for (ki=0; ki<ng[k]; ki++) {
		zerogradGaussian(&gdval);
		for (ni=0; ni<ngp[n]; ni++)
		  for (li=0; li<ng[l]; li++)
		    for (mi=0; mi<ngp[m]; mi++)

		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &Gp[idxp[m]+mi], &Gp[idxp[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idxp[n]+ni)*ngauss],
			     &dGaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &gval, &gdval);

		addmulttogradGaussian(&dval[idx[k]+ki], &gdval, ooa*ovl);
	      }
	      val += 0.5*gval*ooa*ovl;

	      for (s=0; s<A; s++) {
		ooa = o[m+s*A]*o[n+l*A]-o[n+s*A]*o[m+l*A];
		for (si=0; si<ng[s]; si++)
		  addmulttogradGaussian(&dval[idx[s]+si], 
					&dno[(idx[s]+si)+k*ngauss],
					-gval*ooa*ovl);
	      }
	    }
      }	

  grad->val += val;

  // this term is comming from the derivative of the overlap
  for (k=0; k<A; k++)
    for (ki=0; ki<ng[k]; ki++)
      addmulttogradGaussian(&dval[idx[k]+ki], &dno[(idx[k]+ki)+k*ngauss], val);
}


void calcgradSlaterDetTBMErowcol(const SlaterDet* Q, const SlaterDetAux* X,
				 const gradSlaterDetAux* dX,
				 const gradTwoBodyOperator* op, 
				 gradSlaterDet* grad,
				 int k, int l)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* ng = Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double* val = &grad->val; gradGaussian* dval = grad->gradval;

  int m,n,s, ki,li,mi,ni,si;
  complex double gval;
  gradGaussian gdval;
  complex double ooa;

  // zero gradient, as each k,l component is calculated independendly
  *val = 0.0;
  for (ki=0; ki<ngauss; ki++)
    zerogradGaussian(&dval[ki]);

  for (n=0; n<A; n++)
    // for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idx[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  // for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idx[m]*ngauss].T) {

	      ooa = o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A];

	      gval = 0.0;

	      for (ki=0; ki<ng[k]; ki++) {
		zerogradGaussian(&gdval);
		for (ni=0; ni<ng[n]; ni++)
		  for (li=0; li<ng[l]; li++)
		    for (mi=0; mi<ng[m]; mi++)

		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &G[idx[m]+mi], &G[idx[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idx[n]+ni)*ngauss],
			     &dGaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &gval, &gdval);

		addmulttogradGaussian(&dval[idx[k]+ki], &gdval, ooa);
	      }
	      *val += 0.5*gval*ooa;

	      for (s=0; s<A; s++) {
		ooa = o[m+s*A]*o[n+l*A] - o[n+s*A]*o[m+l*A];
		for (si=0; si<ng[s]; si++)
		  addmulttogradGaussian(&dval[idx[s]+si], 
					&dno[(idx[s]+si)+k*ngauss],
					-gval*ooa);
	      }
	    }
      }
}


void calcgradSlaterDetTBMEodrowcol(const SlaterDet* Q, const SlaterDet* Qp, 
				   const SlaterDetAux* X,
				   const gradSlaterDetAux* dX,
				   const gradTwoBodyOperator* op, 
				   gradSlaterDet* grad,
				   int k, int l)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx = Q->idx; int* idxp = Qp->idx; 
  int* ng = Q->ng; int* ngp = Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  gradGaussianAux* dGaux=dX->dGaux;
  gradGaussian* dno=dX->dno;
  complex double val; gradGaussian* dval = grad->gradval;
  complex double ovl = X->ovlap;

  int m,n,s, ki,li,mi,ni,si;
  complex double gval;
  gradGaussian gdval;
  complex double ooa;

  // zero gradient
  grad->val = 0.0;
  for (si=0; si<ngauss; si++)
    zerogradGaussian(&dval[si]);

  val = 0.0;
  for (n=0; n<A; n++)
    // for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idxp[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  // for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idxp[m]*ngauss].T) {

	      ooa = o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A];

	      gval = 0.0;

	      for (ki=0; ki<ng[k]; ki++) {
		zerogradGaussian(&gdval);
		for (ni=0; ni<ngp[n]; ni++)
		  for (li=0; li<ng[l]; li++)
		    for (mi=0; mi<ngp[m]; mi++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &Gp[idxp[m]+mi], &Gp[idxp[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idxp[n]+ni)*ngauss],
			     &dGaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &gval, &gdval);

		addmulttogradGaussian(&dval[idx[k]+ki], &gdval, ooa*ovl);
	      }
	      val += 0.5*gval*ooa*ovl;

	      for (s=0; s<A; s++) {
		ooa = o[m+s*A]*o[n+l*A]-o[n+s*A]*o[m+l*A];
		for (si=0; si<ng[s]; si++)
		  addmulttogradGaussian(&dval[idx[s]+si], 
					&dno[(idx[s]+si)+k*ngauss],
					-gval*ooa*ovl);
	      }
	    }
      }	

  grad->val += val;

  // this term is comming from the derivative of the overlap
  for (s=0; s<A; s++)
    for (si=0; si<ng[s]; si++)
      addmulttogradGaussian(&dval[idx[s]+si], &dno[(idx[s]+si)+s*ngauss], val);
}
