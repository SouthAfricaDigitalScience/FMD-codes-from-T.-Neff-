/**

   \file ProjectedDensityMatrixHO.c

   calculate the projected density matrix in
   harmonic oscillator basis


   (c) 2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "HOBasis.h"
#include "DensityMatrixHO.h"

#include "ProjectedDensityMatrixHO.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MIN(x,y) ((x)<(y) ? (x) : (y))


// make it easy for now: only calculate diagonal


// isoscalar, scalar density matrix
ManyBodyOperator OpDiagonalDensityMatrixHO = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcDiagonalDensityMatrixHOod
};


void initOpDiagonalDensityMatrixHO(DensityMatrixHOPar* par)
{
  char* opname = malloc(40*sizeof(char));
  sprintf(opname, "DiagonalDensityMatrixHO-%d-%05.2f", 
	  par->nmax, hbc*par->omega);
  
  OpDiagonalDensityMatrixHO.name = opname;
  OpDiagonalDensityMatrixHO.dim = (par->nmax+1)*(par->nmax+2);
  OpDiagonalDensityMatrixHO.size = OpDiagonalDensityMatrixHO.dim;
  OpDiagonalDensityMatrixHO.par = par;
}


// workspace

complex double rho[4*HODIMMAX*4*HODIMMAX];


void calcDiagonalDensityMatrixHOod(DensityMatrixHOPar* par,
				   const SlaterDet* Q, const SlaterDet* Qp,
				   const SlaterDetAux* X,
				   complex double* rhodiag)
{
  int dim = dimHOBasis(par->nmax);
  int idx;

  calcDensityMatrixHOod(par, Q, Qp, X, rho);

  int orbit, ixi, xi, N, n, l, twoj, twom; 

  orbit = 0;
  for (N=0; N<=par->nmax; N++) {

    for (ixi=0; ixi<=1; ixi++) {
      xi = ixi ? -1 : +1;

      for (n=N/2; n>=0; n--) {
	l=N-2*n;
	for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2) {
	  rhodiag[orbit] = 0.0;
	  for (twom=-twoj; twom<=twoj; twom=twom+2) {
	    idx = HOxinljmidx(xi, n, l, twoj, twom);
	    rhodiag[orbit] += rho[idx+idx*dim];
	  }
	  orbit++;
	}	
      }
    }
  }
}
 

void writeprojectedDiagonalDensityMatrixHO(FILE* fp,
					   const Projection* P,
					   const DensityMatrixHOPar* dmpar,
					   int j, int pi, int a,
					   void* dmatrix,
					   const Eigenstates* E)
{
  int nmax=dmpar->nmax;
  int dim=(dmpar->nmax+1)*(dmpar->nmax+2);
  int odd=P->odd;
  int jmax=P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int idx = E->index[ipj][a];
  double norm = E->norm[ipj][idx];
  double H = E->v[ipj][idx];

  complex double (**dm)[dim] = dmatrix;
  complex double *noccu = dm[ipj][idx];

  char prefix[8];
  if(odd) sprintf(prefix, "[%d/2%c]", j, pi ? '-' : '+'); 
  else    sprintf(prefix, "[%d%c]", j/2, pi ? '-' : '+'); 


  fprintf(fp, "\n%s (%8.3f MeV) occupation numbers: \n",
	  prefix, hbc*H);

  int i;
  int orbit;
  double nshell, norbit[HONMAX+1];
  int ixi, xi, N, n, l, twoj;

  i=0;
  for (N=0; N<=nmax; N++) {

    fprintf(fp, "\n   N=%d\t\t", N);
    for (n=N/2; n>=0; n--) {
      l=N-2*n;
      for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2)
	fprintf(fp, "%s\t", labelHOorbit(n, l, twoj));
    }
    fprintf(fp, "\n");

    for (ixi=0; ixi<=1; ixi++) {
      xi = ixi ? -1 : +1;
      fprintf(fp, "%c  ", ixi ? 'n' : 'p');

      nshell = 0.0;
      orbit = 0;
      for (n=N/2; n>=0; n--) {
	l=N-2*n;
	for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2) {
	  norbit[orbit] = noccu[i]/norm;
	  i++;
	  nshell += norbit[orbit];
	  orbit++;
	}
      }

      fprintf(fp, "%6.3f\t", nshell);
      for (orbit=0; orbit<=N; orbit++)
	fprintf(fp, "%6.3f\t", norbit[orbit]);
      fprintf(fp, "\n");
    }
  }
  
  double ntot = 0.0;
  for (i=0; i<dim; i++)
    ntot += noccu[i]/norm;

  fprintf(fp, "\nTotal Occupation: %6.3f\n", ntot);

}


#define NMAXOUTPUT 3

void writeprojectedDiagonalDensityMatricesHO(FILE* fp,
					     const Projection* P,
					     const DensityMatrixHOPar* dmpar,
					     void* dmatrix,
					     const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;
  int p,j,i;
  int ngood;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {

      ngood = E->ngood[idxpij(jmax,p,j)];

      for (i=0; i<MIN(ngood, NMAXOUTPUT); i++) {

	writeprojectedDiagonalDensityMatrixHO(fp, P, dmpar, j, p, i, dmatrix, E);
      
      }
    }

}


void writeprojectedDiagonalTransitionDensityMatrixHO(FILE* fp,
						     const Projection* P,
						     const DensityMatrixHOPar* dmpar,
						     int j, int pi, 
						     int af, int ai,
						     void* dmatrix,
						     const Eigenstates* E)
{
  int nmax=dmpar->nmax;
  int dim=(dmpar->nmax+1)*(dmpar->nmax+2);
  int odd=P->odd;
  int jmax=P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int idxf = E->index[ipj][af];
  int idxi = E->index[ipj][ai];
  double normf = E->norm[ipj][idxf];
  double normi = E->norm[ipj][idxi];
  double Hf = E->v[ipj][idxf];
  double Hi = E->v[ipj][idxi];

  complex double (****dm)[dim] = dmatrix;
  complex double *transm = dm[ipj][ipj][idxf][idxi];

  complex double phase;

  // get rid of possible phase factors
  phase = transm[0]/cabs(transm[0]);

  char prefix[8];
  if(odd) sprintf(prefix, "[%d/2%c]", j, pi ? '-' : '+'); 
  else    sprintf(prefix, "[%d%c]", j/2, pi ? '-' : '+'); 


  fprintf(fp, "\n%s (%8.3f MeV -> %8.3f MeV) transition numbers: \n",
	  prefix, hbc*Hi, hbc*Hf);

  int i;
  int orbit;
  double nshell, norbit[HONMAX+1];
  int ixi, xi, N, n, l, twoj;

  i=0;
  for (N=0; N<=nmax; N++) {

    fprintf(fp, "\n   N=%d\t\t", N);
    for (n=N/2; n>=0; n--) {
      l=N-2*n;
      for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2)
	fprintf(fp, "%s\t", labelHOorbit(n, l, twoj));
    }
    fprintf(fp, "\n");

    for (ixi=0; ixi<=1; ixi++) {
      xi = ixi ? -1 : +1;
      fprintf(fp, "%c  ", ixi ? 'n' : 'p');

      nshell = 0.0;
      orbit = 0;
      for (n=N/2; n>=0; n--) {
	l=N-2*n;
	for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2) {
	  norbit[orbit] = creal(transm[i]/phase)/sqrt(normi*normf);
	  i++;
	  nshell += norbit[orbit];
	  orbit++;
	}
      }

      fprintf(fp, "%6.3f\t", nshell);
      for (orbit=0; orbit<=N; orbit++)
	fprintf(fp, "%6.3f\t", norbit[orbit]);
      fprintf(fp, "\n");
    }
  }
  
  double ntot = 0.0;
  for (i=0; i<dim; i++)
    ntot += creal(transm[i]/phase)/sqrt(normi*normf);

  fprintf(fp, "\nTotal Strength: %6.3f\n", ntot);

}
