/**

  \file SlaterDet.c

  Slater determinant build of Gaussians


  (c) 2003 Thomas Neff

*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"
#include "SpatialOrientation.h"
#include "numerics/cmath.h"
#include "numerics/cmat.h"
#include "numerics/rotationmatrices.h"
#include "numerics/lapack.h"

#include "misc/utils.h"

#define SQR(x) (x)*(x)


void allocateSlaterDet(SlaterDet* Q, int A)
{
  Q->A = A;
  Q->idx = malloc(A*sizeof(int));
  Q->ng = malloc(A*sizeof(int));
  Q->G = malloc(MAXNG*A*sizeof(Gaussian));
}


void initSlaterDet(const SlaterDet* Q, SlaterDet* Qp)
{
  assert(Q->ngauss <= MAXNG*Q->A);
  Qp->A = Q->A;
  Qp->ngauss = Q->ngauss;

  Qp->idx = malloc(Qp->A*sizeof(int));
  Qp->ng = malloc(Qp->A*sizeof(int));
  Qp->G = malloc(MAXNG*Qp->A*sizeof(Gaussian));
}


void freeSlaterDet(SlaterDet* Q)
{
  free(Q->idx);
  free(Q->ng);
  free(Q->G);
}


void writeSlaterDet(FILE* fp, const SlaterDet* Q)
{
  int k, i, idx;

  fprintf(fp, "<SlaterDet>\n");
  fprintf(fp, "%d %d %d\n", Q->A, Q->Z, Q->N);
  fprintf(fp, "%d\n", Q->ngauss);
  for (k=0; k<Q->A; k++)
    for (i=0; i<Q->ng[k]; i++) {
      idx = Q->idx[k]+i;
      fprintf(fp,
	      "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	    k, Q->G[idx].xi, 
	    creal(Q->G[idx].chi[0]), cimag(Q->G[idx].chi[0]), 
	    creal(Q->G[idx].chi[1]), cimag(Q->G[idx].chi[1]),
	    creal(Q->G[idx].a), cimag(Q->G[idx].a),
	    creal(Q->G[idx].b[0]), cimag(Q->G[idx].b[0]),
	    creal(Q->G[idx].b[1]), cimag(Q->G[idx].b[1]),
	    creal(Q->G[idx].b[2]), cimag(Q->G[idx].b[2]));
  }

  fprintf(fp, "</SlaterDet>\n");
}


#define BUFSIZE 255

int readSlaterDet(FILE* fp, SlaterDet* Q)
{
  char buf[BUFSIZE];
  int i;
  
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<SlaterDet>", 11) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <SlaterDet>\n");
    return -1;
  }

  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d %d %d", &Q->A, &Q->Z, &Q->N) != 3) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }
  if (Q->A != Q->Z+Q->N) {
    fprintf(stderr, "proton and neutron number don't sum up in line\n>>%s<<\n",
	    stripstr(buf, "\n"));
    return -1;
  }
  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d", &Q->ngauss) != 1) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }

  Q->idx = malloc(Q->A*sizeof(int));
  Q->ng = malloc(Q->A*sizeof(int));
  for (i=0; i<Q->A; i++) {
    Q->idx[i] = 0;
    Q->ng[i] = 0;
  }
  Q->G = malloc(MAXNG*Q->A*sizeof(Gaussian));

  int k, xi;
  double chi0re, chi0im, chi1re, chi1im;
  double are, aim, b0re, b0im, b1re, b1im, b2re, b2im; 
  for (i=0; i<Q->ngauss; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "%d %d "
	       "(%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf)", 
	       &k, &xi, &chi0re, &chi0im, &chi1re, &chi1im, 
	       &are, &aim, &b0re, &b0im, &b1re, &b1im, &b2re, &b2im) != 14) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    Q->ng[k]++;
    Q->G[i].xi = xi;
    Q->G[i].chi[0] = chi0re + I*chi0im;
    Q->G[i].chi[1] = chi1re + I*chi1im;
    Q->G[i].a = are + I*aim;
    Q->G[i].b[0] = b0re + I*b0im;
    Q->G[i].b[1] = b1re + I*b1im;
    Q->G[i].b[2] = b2re + I*b2im;
  }	
  for (k=1; k<Q->A; k++) 
    Q->idx[k] = Q->idx[k-1]+Q->ng[k-1];
	 
  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</SlaterDet>", 12)) {
    fprintf(stderr, "did't find </SlaterDet>\n");
    return -1;
  }

  return 0;
}
  

/// read Slater determinant Q from file
/// Q will be initialized
int readSlaterDetfromFile(SlaterDet* Q, const char* fname)
{
  FILE* fp;

  if (!(fp = fopen(fname, "r"))) {  
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    return -1;
  }

  fprintf(stderr, "... reading SlaterDet from file %s\n", fname);

  int err = readSlaterDet(fp, Q);
  fclose(fp);

  return err;
}


void copySlaterDet(const SlaterDet* Q, SlaterDet* Qp) 
{
  int i;
  
  assert(Qp->A >= Q->A);

  Qp->A = Q->A;
  Qp->Z = Q->Z; Qp->N = Q->N;
  Qp->ngauss = Q->ngauss;
  for (i=0; i<Q->A; i++) {
    Qp->idx[i] = Q->idx[i];
    Qp->ng[i] = Q->ng[i];
  }
  for (i=0; i<Q->ngauss; i++)
    Qp->G[i] = Q->G[i];
}


void cloneSlaterDet(const SlaterDet* Q, SlaterDet* Qp)
{
  Qp->A = Q->A;
  Qp->ngauss = Q->ngauss;

  Qp->idx = malloc(Qp->A*sizeof(int));
  Qp->ng = malloc(Qp->A*sizeof(int));
  Qp->G = malloc(MAXNG*Qp->A*sizeof(Gaussian));

  copySlaterDet(Q, Qp);
}


void joinSlaterDets(const SlaterDet* Qa, const SlaterDet* Qb, SlaterDet* Q)
{
  Q->A = Qa->A+Qb->A;
  Q->Z = Qa->Z+Qb->Z;
  Q->N = Qa->N+Qb->N;
  Q->ngauss = Qa->ngauss+Qb->ngauss;

  Q->idx = malloc(Q->A*sizeof(int));
  Q->ng = malloc(Q->A*sizeof(int));
  Q->G = malloc(MAXNG*Q->A*sizeof(Gaussian));

  int i;
  // copy Qa
  for (i=0; i<Qa->A; i++) {
    Q->idx[i] = Qa->idx[i];
    Q->ng[i] = Qa->ng[i];
  }
  for (i=0; i<Qa->ngauss; i++)
    Q->G[i] = Qa->G[i];

  // copy Qb
  for (i=0; i<Qb->A; i++) {
    Q->idx[i+Qa->A] = Qb->idx[i]+Qa->ngauss;
    Q->ng[i+Qa->A] = Qb->ng[i];
  }
  for (i=0; i<Qb->ngauss; i++) 
    Q->G[i+Qa->ngauss] = Qb->G[i];
}


void normalizeSlaterDet(SlaterDet* Q, SlaterDetAux* X)
{
  calcSlaterDetAuxod(Q, Q, X);
  double norm = sqrt(creal(X->ovlap));
  double normsp = pow(norm, 1.0/Q->A);

  int i;
  for (i=0; i<Q->ngauss; i++) {
    Q->G[i].chi[0] /= normsp;
    Q->G[i].chi[1] /= normsp;
  }
}


void moveSlaterDet(SlaterDet* Q, double d[3])
{
  int i;

  for (i=0; i<Q->ngauss; i++)
    moveGaussian(&Q->G[i], d);
}


void boostSlaterDet(SlaterDet* Q, double v[3])
{
  int i;

  for (i=0; i<Q->ngauss; i++)
    boostGaussian(&Q->G[i], v);
}


void moveboostSlaterDet(SlaterDet* Q, SlaterDetAux* X)
{
  int i;

  // move SlaterDet to origin in phase space

  double Xcm[3], Vcm[3];
  calcSlaterDetAux(Q, X);
  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);  

  for (i=0; i<3; i++)
    Vcm[i] *= -1;
  for (i=0; i<3; i++)
    Xcm[i] *= -1;

  boostSlaterDet(Q, Vcm);
  moveSlaterDet(Q, Xcm);
}


void rotateSlaterDet(SlaterDet* Q, double a, double b, double g)
{
  double euler[3] = { a, b, g };

  complex double R2[2][2];
  rotatemat2(euler, R2);

  double R3[3][3];
  rotatemat3(euler, R3);

  int i;
  for (i=0; i<Q->ngauss; i++)
    rotateGaussian(&Q->G[i], R3, R2);
}


void moveboostorientSlaterDet(SlaterDet* Q, SlaterDetAux* X)
{
  int i;

  // move SlaterDet to origin in phase space
  // orient so that min moment of inertia around the z-axes

  double Xcm[3], Vcm[3];
  calcSlaterDetAux(Q, X);
  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);  

  for (i=0; i<3; i++)
    Vcm[i] *= -1;
  for (i=0; i<3; i++)
    Xcm[i] *= -1;

  boostSlaterDet(Q, Vcm);
  moveSlaterDet(Q, Xcm);

  // orient the nucleus along the main axes

  double alpha, beta, gamma;
  calcSlaterDetAux(Q, X);
  calcOrientedOrientation(Q, X, &alpha, &beta, &gamma);
  rotateSlaterDet(Q, -gamma, -beta, -alpha);
}


void orientSlaterDet(SlaterDet* Q, SlaterDetAux* X)
{
  // orient the nucleus along the main axes

  double alpha, beta, gamma;
  calcSlaterDetAux(Q, X);
  calcOrientedOrientation(Q, X, &alpha, &beta, &gamma);
  rotateSlaterDet(Q, -gamma, -beta, -alpha);
}


void invertSlaterDet(SlaterDet* Q)
{
  int i;

  for (i=0; i<Q->ngauss; i++)
    invertGaussian(&Q->G[i]);
}


void invertcmSlaterDet(SlaterDet* Q, SlaterDetAux* X)
{
  double Xcm[3], Vcm[3];
  int i;

  calcSlaterDetAux(Q, X);
  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  for (i=0; i<3; i++) {
    Xcm[i] *= -1;
    Vcm[i] *= -1;
  }
  boostSlaterDet(Q, Vcm);
  moveSlaterDet(Q, Xcm);

  invertSlaterDet(Q);

  for (i=0; i<3; i++) {
    Xcm[i] *= -1;
    Vcm[i] *= -1;
  }
  moveSlaterDet(Q, Xcm);
  boostSlaterDet(Q, Vcm);

}


void timerevertSlaterDet(SlaterDet* Q)
{
  int i;

  for (i=0; i<Q->ngauss; i++)
    timerevertGaussian(&Q->G[i]);
}


void scaleSlaterDet(SlaterDet* Q, double kappa)
{
  int i;

  for (i=0; i<Q->ngauss; i++)
    scaleGaussian(&Q->G[i], kappa);
}


void allocateSlaterDetAux(SlaterDetAux* X, int A)
{
  X->Gaux = malloc(SQR(MAXNG*A)*sizeof(GaussianAux));
  X->n = initcmat(A);
  X->o = initcmat(A);
}


// space is reserved for up to MAXNG Gaussians per nucleon
void initSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X)
{
  assert(Q->ngauss <= MAXNG*Q->A);
  X->Gaux = malloc(SQR(MAXNG*Q->A)*sizeof(GaussianAux));
  X->n = initcmat(Q->A);
  X->o = initcmat(Q->A);
}  
  

void freeSlaterDetAux(SlaterDetAux* X)
{
  free(X->Gaux);
  free(X->n);
  free(X->o);
}


// hermiticity of matrices is not exploited
void calcSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;

  int k,l,ki,li;

  for (l=0; l<A; l++)
    for (li=0; li<ng[l]; li++)
      for (k=0; k<A; k++)
	for (ki=0; ki<ng[k]; ki++)
	  calcGaussianAux(&G[idx[k]+ki], &G[idx[l]+li], 
			  &Gaux[(idx[k]+ki)+(idx[l]+li)*ngauss]);


  for (l=0; l<A; l++) 
    for (k=0; k<A; k++) {
      X->n[k+l*A] = 0.0;
      for (li=0; li<ng[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  X->n[k+l*A] += Gaux[(idx[k]+ki)+(idx[l]+li)*ngauss].Q;
    }

  char UPLO='U';
  int info;

  copycmat(A, X->n, X->o);
  FORTRAN(zpotrf)(&UPLO, &A, X->o, &A, &info);
  FORTRAN(zpotri)(&UPLO, &A, X->o, &A, &info);
  mirrorcmat(A, X->o);
}


///
void calcSlaterDetAuxod(const SlaterDet* Q, const SlaterDet* Qp,
			SlaterDetAux* X)
{
  assert(Q->A == Qp->A && Q->Z == Qp->Z && Q->N == Qp->N);

  int A=Q->A; int ngauss=Q->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* n=X->n;

  int k,l,ki,li;

  for (l=0; l<A; l++)
    for (li=0; li<ngp[l]; li++)
      for (k=0; k<A; k++)
	for (ki=0; ki<ng[k]; ki++)
	  calcGaussianAux(&G[idx[k]+ki], &Gp[idxp[l]+li], 
			  &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss]);


  for (l=0; l<A; l++) 
    for (k=0; k<A; k++) {
      n[k+l*A] = 0.0;
      for (li=0; li<ngp[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  n[k+l*A] += Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss].Q;
    }


  int ipiv[A];
  int lwork=A*A;
  complex double work[lwork];
  int info;
  
  copycmat(A, X->n, X->o);
  FORTRAN(zgetrf)(&A, &A, X->o, &A, ipiv, &info);
  FORTRAN(zdet)(X->o, &A, &A, ipiv, &X->ovlap);
  FORTRAN(zgetri)(&A, X->o, &A, ipiv, work, &lwork, &info);

  if (info)
    fprintf(stderr, "calcSlaterDetAuxod: overlap matrix singular !\n");
}


// sort of a hack
// calculate cofactors by using the svd 
// assume rank A-1 for overlap matrix
// ovlap is set to 1 
// only use for off-diagonal matrix elements of One-Body Operators
// broken for Two-body Operators
void calcSlaterDetAuxodsingular(const SlaterDet* Q, const SlaterDet* Qp,
				SlaterDetAux* X)
{
  int A=Q->A; int ngauss=Q->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* n=X->n;
  complex double* o=X->o;

  int k,l,ki,li;

  for (l=0; l<A; l++)
    for (li=0; li<ngp[l]; li++)
      for (k=0; k<A; k++)
	for (ki=0; ki<ng[k]; ki++)
	  calcGaussianAux(&G[idx[k]+ki], &Gp[idxp[l]+li], 
			  &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss]);


  for (l=0; l<A; l++) 
    for (k=0; k<A; k++) {
      n[k+l*A] = 0.0;
      for (li=0; li<ngp[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  n[k+l*A] += Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss].Q;
    }

  // calculate SVD of overlap matrix
  complex double U[A*A];
  complex double Vh[A*A];
  double S[A];
  {
    const char jobu='A', jobvt='A';
    const int lwork=5*A;
    complex double work[lwork];
    double rwork[lwork];
    int info;

    FORTRAN(zgesvd)(&jobu, &jobvt, &A, &A, n, &A, S,
		    U, &A, Vh, &A, work, &lwork, rwork, &info);
  }

  // calculate phase detVhU using LU decomposition
  complex double detVhU;
  {
    complex double VhU[A*A];
    int ipiv[A];
    int info;

    multcmat(Vh, U, VhU, A);
    FORTRAN(zgetrf)(&A, &A, VhU, &A, ipiv, &info);
    FORTRAN(zdet)(VhU, &A, &A, ipiv, &detVhU);
  }

  int i,j;
  double prodS;

  prodS = 1.0;
  for (j=0; j<A-1; j++)
    prodS *= S[j];
  
  for (l=0; l<A; l++)
    for (k=0; k<A; k++)
      o[k+l*A] = conj(Vh[A-1+k*A])*prodS*conj(U[l+(A-1)*A])*detVhU;

  X->ovlap = 1.0;
}


void calcSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
		       const OneBodyOperator* op, double val[])
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;

  int k,l,ki,li;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (l=0; l<A; l++) {
    for (k=0; k<l; k++) 
      if (!op->opt || Gaux[idx[k]+idx[l]*ngauss].T) {
      for (i=0; i<op->dim; i++) 
	gval[i] = 0.0;

      for (li=0; li<ng[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  op->me(op->par, 
		 &G[idx[k]+ki], &G[idx[l]+li], 
		 &Gaux[(idx[k]+ki)+(idx[l]+li)*ngauss], 
		 gval);

      for (i=0; i<op->dim; i++)
	val[i] += 2.0*gval[i]*o[l+k*A];
    }	
    
    for (i=0; i<op->dim; i++) 
      gval[i] = 0.0;
    for (li=0; li<ng[l]; li++)
      for (ki=0; ki<ng[l]; ki++)
	op->me(op->par,
	       &G[idx[l]+ki], &G[idx[l]+li], 
	       &Gaux[(idx[l]+ki)+(idx[l]+li)*ngauss], 
	       gval);	 

    for (i=0; i<op->dim; i++)
      val[i] += creal(gval[i])*creal(o[k+l*A]);
  }	

  free(gval);
}


void calcSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
		       const TwoBodyOperator* op, double val[])
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;

  int k,l,m,n, ki,li,mi,ni;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (n=0; n<A; n++)
    for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idx[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  for (k=0; k<l; k++)
	    if (!op->opt || Gaux[idx[k]+idx[m]*ngauss].T) {

	      for (i=0; i<op->dim; i++) 
		gval[i] = 0.0;

	      for (ni=0; ni<ng[n]; ni++)
		for (li=0; li<ng[l]; li++)
		  for (mi=0; mi<ng[m]; mi++)
		    for (ki=0; ki<ng[k]; ki++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &G[idx[m]+mi], &G[idx[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idx[n]+ni)*ngauss],
			     gval);

	      for (i=0; i<op->dim; i++)
		val[i] += gval[i]*
		  (o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A]);
	    }	
      }
  
  free(gval);
}


void calcSlaterDetTBMErowcol(const SlaterDet* Q, const SlaterDetAux* X,
			     const TwoBodyOperator* op, double val[], 
			     int k, int l)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;

  int m,n, ki,li,mi,ni;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (n=0; n<A; n++)
    // for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idx[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  // for (k=0; k<l; k++)
	    if (!op->opt || Gaux[idx[k]+idx[m]*ngauss].T) {

	      for (i=0; i<op->dim; i++) 
		gval[i] = 0.0;

	      for (ni=0; ni<ng[n]; ni++)
		for (li=0; li<ng[l]; li++)
		  for (mi=0; mi<ng[m]; mi++)
		    for (ki=0; ki<ng[k]; ki++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &G[idx[m]+mi], &G[idx[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idx[n]+ni)*ngauss],
			     gval);

	      for (i=0; i<op->dim; i++)
		val[i] += gval[i]*
		  (o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A]);
	    }	
      }	

  free(gval);
}


void calcSlaterDetOBMEod(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X,
			 const OneBodyOperator* op, complex double val[])
{
  int A=Q->A; int ngauss=Q->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  complex double ovl = X->ovlap;

  int k,l, ki,li;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (l=0; l<A; l++)
    for (k=0; k<A; k++) 
      if (!op->opt || Gaux[idx[k]+idxp[l]*ngauss].T) {

	for (i=0; i<op->dim; i++) 
	  gval[i] = 0.0;

	for (li=0; li<ngp[l]; li++)
	  for (ki=0; ki<ng[k]; ki++)
	    op->me(op->par, 
		   &G[idx[k]+ki], &Gp[idxp[l]+li], 
		   &Gaux[(idx[k]+ki)+(idxp[l]+li)*ngauss], 
		   gval);

	for (i=0; i<op->dim; i++)
	  val[i] += gval[i]*o[l+k*A]*ovl;
      }		
  
  free(gval);
}	


void calcSlaterDetTBMEod(const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X,
			 const TwoBodyOperator* op, complex double val[])
{
  int A=Q->A; int ngauss=Q->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  complex double ovl=X->ovlap;

  int k,l,m,n, ki,li,mi,ni;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (n=0; n<A; n++)
    for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idxp[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idxp[m]*ngauss].T) {

	      for (i=0; i<op->dim; i++) 
		gval[i] = 0.0;

	      for (ni=0; ni<ngp[n]; ni++)
		for (li=0; li<ng[l]; li++)
		  for (mi=0; mi<ngp[m]; mi++)
		    for (ki=0; ki<ng[k]; ki++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &Gp[idxp[m]+mi], &Gp[idxp[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idxp[n]+ni)*ngauss],
			     gval);

	      for (i=0; i<op->dim; i++)
		val[i] += 0.5*gval[i]*
		  (o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A])*ovl;
	    }	
      }

  free(gval);
}


void calcSlaterDetTBMEodrowcol(const SlaterDet* Q, const SlaterDet* Qp,
			       const SlaterDetAux* X,
			       const TwoBodyOperator* op, complex double val[],
			       int k, int l)
{
  int A=Q->A; int ngauss=Q->ngauss;
  int* idx=Q->idx; int* idxp=Qp->idx; 
  int* ng=Q->ng; int* ngp=Qp->ng;
  Gaussian* G=Q->G; Gaussian* Gp=Qp->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  complex double ovl=X->ovlap;

  int m,n, ki,li,mi,ni;
  int i;
  complex double *gval = malloc(op->dim*sizeof(complex double));

  for (i=0; i<op->dim; i++)
    val[i] = 0.0;

  for (n=0; n<A; n++)
    // for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idxp[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  // for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idxp[m]*ngauss].T) {

	      for (i=0; i<op->dim; i++) 
		gval[i] = 0.0;

	      for (ni=0; ni<ngp[n]; ni++)
		for (li=0; li<ng[l]; li++)
		  for (mi=0; mi<ngp[m]; mi++)
		    for (ki=0; ki<ng[k]; ki++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &Gp[idxp[m]+mi], &Gp[idxp[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idxp[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idxp[n]+ni)*ngauss],
			     gval);

	      for (i=0; i<op->dim; i++)
		val[i] += gval[i]*
		  (o[m+k*A]*o[n+l*A]-o[n+k*A]*o[m+l*A])*ovl;
	    }	
      }

  free(gval);
}



void calcSlaterDetOBHFMEs(const SlaterDet* Q, const SlaterDetAux* X,
			  const OneBodyOperator* op, void* val)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double (*mes)[op->dim] = val;

  int k,l,ki,li;
  int i;
  complex double gval[op->dim];

  for (l=0; l<A; l++)
    for (k=0; k<A; k++)
      for (i=0; i<op->dim; i++)
	mes[k+l*A][i] = 0.0;

  for (l=0; l<A; l++) {
    for (k=0; k<A; k++) 
      if (!op->opt || Gaux[idx[k]+idx[l]*ngauss].T) {
      for (i=0; i<op->dim; i++) 
	gval[i] = 0.0;

      for (li=0; li<ng[l]; li++)
	for (ki=0; ki<ng[k]; ki++)
	  op->me(op->par, 
		 &G[idx[k]+ki], &G[idx[l]+li], 
		 &Gaux[(idx[k]+ki)+(idx[l]+li)*ngauss], 
		 gval);

      for (i=0; i<op->dim; i++)
	mes[k+l*A][i] += gval[i];
    }	
  }	
}


void calcSlaterDetTBHFMEs(const SlaterDet* Q, const SlaterDetAux* X,
			  const TwoBodyOperator* op, void* val)
{
  int A=Q->A; int ngauss=Q->ngauss; 
  int* idx=Q->idx; int* ng=Q->ng;
  Gaussian* G=Q->G;
  GaussianAux* Gaux=X->Gaux;
  complex double* o=X->o;
  complex double (*mes)[op->dim] = val;

  int k,l,m,n, ki,li,mi,ni;
  int i;
  complex double gval[op->dim];

  for (m=0; m<A; m++)
    for (k=0; k<A; k++)
      for (i=0; i<op->dim; i++)
	mes[k+m*A][i] = 0.0;

  for (n=0; n<A; n++)
    for (l=0; l<A; l++)
      if (!op->opt || Gaux[idx[l]+idx[n]*ngauss].T) {
	for (m=0; m<A; m++)
	  for (k=0; k<A; k++)
	    if (!op->opt || Gaux[idx[k]+idx[m]*ngauss].T) {

	      for (i=0; i<op->dim; i++) 
		gval[i] = 0.0;

	      for (ni=0; ni<ng[n]; ni++)
		for (li=0; li<ng[l]; li++)
		  for (mi=0; mi<ng[m]; mi++)
		    for (ki=0; ki<ng[k]; ki++)
		      op->me(op->par, 
			     &G[idx[k]+ki], &G[idx[l]+li], 
			     &G[idx[m]+mi], &G[idx[n]+ni], 
			     &Gaux[(idx[k]+ki)+(idx[m]+mi)*ngauss],
			     &Gaux[(idx[l]+li)+(idx[n]+ni)*ngauss],
			     gval);

	      for (i=0; i<op->dim; i++) {
		mes[k+m*A][i] += gval[i]*o[n+l*A];
		mes[k+n*A][i] -= gval[i]*o[m+l*A];
	      }
	    }	
      }	
}
