/*

  cmat.c

  complex matrices


  (c) 2003,2004 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "cmat.h"
#include "lapack.h"


complex double* initcmat(int n)
{
  int i;
  complex double* A;

  A = (complex double*) malloc(n*n*sizeof(complex double));
  for (i=0; i<n*n; i++)
    A[i] = 0.0;

  return A;
}


void copycmat(int n, const complex double* A, complex double* B)
{
  int i;

  for (i=0; i<n*n; i++)
    B[i] = A[i];
}


void mirrorcmat(int n, complex double* A)
{
  int i,j;

  for (j=0; j<n; j++)
    for (i=0; i<j; i++)
      A[j+i*n] = conj(A[i+j*n]);
}


void hermitizecmat(int n, complex double* A)
{
  int i,j;
  complex double alower, aupper;

  for (i=0; i<n; i++)
    for (j=0; j<=i; j++) {
      aupper = A[i+n*j];
      alower = A[j+n*i];
      A[i+n*j] = 0.5*(aupper+conj(alower));
      A[j+n*i] = 0.5*(alower+conj(aupper));
    }
}


void multcmat(const complex double* A, const complex double* B,
	      complex double* C, int n)
{
  int k,l,m;

  for (l=0; l<n; l++)
    for (k=0; k<n; k++) {
      C[l*n+k] = 0.0;
      for (m=0; m<n; m++)
	C[l*n+k] += A[m*n+k]*B[l*n+m];
    }
}


#define BUFSIZE 4096
void freadivec(FILE* fp, int n, int* v)
{
  int k;
  char buf[BUFSIZE];
  char* c;

  fgets(buf, BUFSIZE, fp);
  
  c=strtok(buf, " ");
  for (k=0; k<n; k++) {
    v[k]=atoi(c);
    c=strtok(NULL, " ");
  }
}


void fprintcvec(FILE* fp, int n, const complex double* V)
{
  int i;

  for (i=0; i<n; i++)
    fprintf(fp, "(%15.8e,%15.8e) ", creal(V[i]), cimag(V[i]));
  fprintf(fp, "\n");
}


void fwritecvecbin(FILE* fp, int n, const complex double* V)
{
  fwrite(V, sizeof (complex double), n, fp);
}
    

#define BUFSIZE 655536
void freadcvec(FILE* fp, int n, complex double* a)
{
  int l;
  char buf[BUFSIZE];
  char *c;
  double ref, imf;

  fgets(buf, BUFSIZE, fp);

  c = strtok(buf, " ,()");
  for (l=0; l<n; l++) {
    ref=atof(c); c=strtok(NULL, " ,()");
    imf=atof(c); c=strtok(NULL, " ,()");
    a[l] = ref + I*imf;
  }
}


void fprintcmat(FILE* fp, int n, const complex double* A)
{
  int i,j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      fprintf(fp, "(%15.8e,%15.8e) ", creal(A[i+j*n]), cimag(A[i+j*n]));
    fprintf(fp, "\n");
  }
}


void fwritecmatbin(FILE* fp, int n, const complex double* A)
{
  fwrite(A, sizeof(complex double), n*n, fp);
}


#define BUFSIZE 655536
void freadcmat(FILE* fp, int n, complex double* A)
{
  int k, l;
  char buf[BUFSIZE];
  char *c;
  double ref, imf;

  for (k=0; k<n; k++) {
    fgets(buf, BUFSIZE, fp);
    c = strtok(buf, " ,()");
    for (l=0; l<n; l++) {
      ref=atof(c); c=strtok(NULL, " ,()");
      imf=atof(c); c=strtok(NULL, " ,()");
      A[k+l*n] = ref + I*imf;
    }
  }
}


// matrix is transposed !
void fprintcmatcols(FILE* fp, int n, int cols, const complex double* A)
{
  int i,j;

  for (i=0; i<cols; i++) {
    for (j=0; j<n; j++)
      fprintf(fp, "(%15.8e,%15.8e) ", creal(A[j+i*n]), cimag(A[j+i*n]));
    fprintf(fp, "\n");
  }
}


// matrix is transposed !
#define BUFSIZE 655536
void freadcmatcols(FILE* fp, int n, int cols, complex double* A)
{
  int k, l;
  char buf[BUFSIZE];
  char *c;
  double ref, imf;

  for (k=0; k<cols; k++) {
    fgets(buf, BUFSIZE, fp);
    c = strtok(buf, " ,()");
    for (l=0; l<n; l++) {
      ref=atof(c); c=strtok(NULL, " ,()");
      imf=atof(c); c=strtok(NULL, " ,()");
      A[l+k*n] = ref + I*imf;
    }
  }
}


void pseudoinverse(const complex double* A, complex double* B,
		   int n, double thresh)
{
  const char jobu='O', jobvt='A';

  complex double *U = malloc(n*n*sizeof(complex double));
  complex double *Vh = malloc(n*n*sizeof(complex double));
  double* S = malloc(n*sizeof(double));

  const int lwork=5*n;
  complex double* work = malloc(lwork*sizeof(complex double));
  double* rwork = malloc(lwork*sizeof(double));

  int info;
  
  copycmat(n, A, U);
  FORTRAN(zgesvd)(&jobu, &jobvt, &n, &n, U, &n, S, 
		  NULL, &n, Vh, &n, work, &lwork, rwork, &info);

  int k,l,m;

  for (l=0; l<n; l++)
    for (k=0; k<n; k++) {
      B[k+l*n] = 0.0;
      m=0; 
      while (m<n && S[m]/S[0] > thresh) {
	B[k+l*n] += conj(Vh[m+k*n])*1/S[m]*conj(U[l+m*n]);
	m++;
      }	
    }

  free(Vh); free(U); free(S);
  free(work); free(rwork);
}


void eigensystem(complex double* A,
		 complex double* a, complex double* V, int n)
{
  const char jobvl='N', jobvr='V';

  const int lwork=5*n;
  complex double* work = malloc(lwork*sizeof(complex double));
  double* rwork = malloc(2*n*sizeof(double));
  int info;
  
  FORTRAN(zgeev)(&jobvl, &jobvr, &n, A, &n, a, NULL, &n, V, &n,
		 work, &lwork, rwork, &info);
  if (info) {
    fprintf(stderr, "eigensystem: zgeev returned with error code: %d\n", info);
  }

  free(work); free(rwork);
}


void generalizedeigensystem(complex double* A, complex double* N, int n,
			    double thresh,
			    complex double* v, complex double* V, int* dim)
{
  const char jobu='O', jobvt='A';

  complex double* U = malloc(n*n*sizeof(complex double));
  complex double* Vh = malloc(n*n*sizeof(complex double));
  double* S = malloc(n*sizeof(double));

  const int lwork=5*n;
  complex double* work = malloc(lwork*sizeof(complex double));
  double* rwork = malloc(lwork*sizeof(double));

  int info;
  
  // perform SVD of overlap matrix
  copycmat(n, N, U);
  FORTRAN(zgesvd)(&jobu, &jobvt, &n, &n, U, &n, S, 
		  NULL, &n, Vh, &n, work, &lwork, rwork, &info);

  // for zero matrices
  if (S[0] == 0.0) {
    *dim = 0;
    return;
  }

  // dimension of subspace
  int m=0;

  while (m<n && S[m]/S[0] > thresh)
    m++;

  // if our subspace is empty nothing to do any more

  if (!m) {
    *dim = 0;
    return;
  }

  // project into subspace
  complex double *AV = malloc(m*n*sizeof(complex double));
  complex double *invNA = malloc(m*m*sizeof(complex double));
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
  complex double *vb = malloc(m*sizeof(complex double));
  complex double *Vb = malloc(m*m*sizeof(complex double));

  FORTRAN(zgeev)(&jobvl, &jobvr, &m, invNA, &m, vb, NULL, &m, Vb, &m,
		 work, &lwork, rwork, &info);

  // imbed solution in full space

  for (i=0; i<m; i++)
    v[i] = vb[i];

  for (k=0; k<m; k++)
    for (i=0; i<n; i++) {
      V[i+k*n] = 0.0;
      for (l=0; l<m; l++)
	V[i+k*n] += conj(Vh[l+i*n])*Vb[l+k*m];
    }

  // // normalize eigenvectors with respect to overlap matrix N
  //
  // double *norm2 = malloc(m*sizeof(double));
  //
  // for (i=0; i<m; i++) {
  //   norm2[i] = 0.0;
  //
  //   for (k=0; k<n; k++)
  //     for (l=0; l<n; l++)
  //       norm2[i] += conj(V[l+i*n])*N[l+k*n]*V[k+i*n];
  // }
  //
  // for (i=0; i<m; i++) 
  //   for (k=0; k<n; k++)
  //     V[k+i*n] /= sqrt(norm2[i]);
  //
  // free(norm2);

  *dim = m;

  free(Vb); free(vb);
  free(AV); free(invNA);
  free(work); free(rwork);
  free(U); free(Vh); free(S);
}


// sort first dim eigenstates, starting with smalles eigenvalue

static int cmpmerit(void* ap, void* bp) 
{
  double a=*(double*) ap;
  double b=*(double*) bp;
  
  return (a<=b ? (a<b ? 1 : 0) : -1);
}

void sorteigenstates(int n, complex double* v, complex double* V, int dim)
{
  int i,k;

  struct merit {
    double val;
    int idx;
  } merits[dim];

  for (i=0; i<dim; i++) {
    merits[i].idx = i;
    merits[i].val = -creal(v[i]);
  }

  qsort(merits, dim, sizeof(struct merit), cmpmerit);

  complex double vsorted[dim];
  complex double Vsorted[n*dim];

  for (i=0; i<dim; i++) {
    vsorted[i] = v[merits[i].idx];
    for (k=0; k<n; k++) 
      Vsorted[k+i*n] = V[k+merits[i].idx*n];
  }

  // copy into original v and V

  for (i=0; i<dim; i++) {
    v[i] = vsorted[i];
    for (k=0; k<n; k++)
      V[k+i*n] = Vsorted[k+i*n];
  }
}
