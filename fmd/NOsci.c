/**

  \file NOsci.c 

  Harmonic Oscillator quanta

  (c) 2007 Thomas Neff

*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "NOsci.h"

#include "numerics/cmath.h"
#include "misc/utils.h"
#include "misc/physics.h"


// NOsci
ManyBodyOperator OpNOsci = {
  name : "NOsci",
  rank : 0,
  pi : 0,
  dim : 2,
  size : 2,
  par : NULL,
  me : calcNOsciod
};


// 0hbw oscillator quanta

int nq[] = {  0,                                 // no nucleons
              0, 0,                              // s-shell
              1, 2, 3, 4, 5, 6,                  // p-shell
              8,10,12,14,16,18,20,22,24,26,      // sd-shell
              29,32,35,38,41,44,47,50,53,56,     // pf-shell
              59,62,65,68,71,74,77,80,83,86 };

static inline double nquanta0(const SlaterDet* Q)
{
  return (1.5*Q->A + nq[Q->Z] + nq[Q->N]);
}

// calculate x^2 and p^2 matrix elements
static void ob_nosci(void* par,
                     const Gaussian*G1, const Gaussian* G2, 
                     const GaussianAux* X, 
                     complex double nosci[2])
{
  // assume Slater determinants are in origin 
  double Xcm[3] = {0, 0, 0};
  double Vcm[3] = {0, 0, 0};

  int i;
  complex double rho2cm, pi2cm;

  rho2cm = 3*X->alpha;
  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  pi2cm = 3*X->lambda;
  for (i=0; i<3; i++)
    pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);
  
  nosci[0] = rho2cm* X->Q;
  nosci[1] = pi2cm* X->Q;
}


void calcNOsciod(void* Par,
                 const SlaterDet* Q, const SlaterDet* Qp,
                 const SlaterDetAux* X, 
                 complex double nosci[2])
{
  OneBodyOperator op_ob_nosci = {dim: 2, opt: 1, par: NULL, me: ob_nosci};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_nosci, nosci);
}


void writeprojectedNOsci(FILE* fp,
                         const Projection* P,
                         const complex double (**nosci)[2],
                         const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i;

  char prefix[8];

  int* idx;
  int ngood;
  complex double* norm;
  complex double *H, (*NO)[2];
  double x2, p2;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

	if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
	else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
	norm = E->norm[idxpij(jmax,p,j)];
	H = E->v[idxpij(jmax,p,j)];
 
	NO = nosci[idxpij(jmax,p,j)];

	fprintf(fp, "\n%s  N     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
	fprintf(fp, "\n%s  H     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
	fprintf(fp, "\n%s  NOsci = ", prefix); 
	for (i=0; i<ngood; i++) {
          x2 = creal(NO[idx[i]][0]/norm[idx[i]]);
          p2 = creal(NO[idx[i]][1]/norm[idx[i]]);
	  fprintf(fp, "   %8.3f", sqrt(x2*p2));
        }
	fprintf(fp, "\n");	  
      }
    }
}
