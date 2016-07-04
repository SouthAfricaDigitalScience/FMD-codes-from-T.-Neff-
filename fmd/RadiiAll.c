/**

  \file RadiiAll.c

  calculate radii


  (c) 2009 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "RadiiAll.h"

#include "misc/physics.h"
#include "numerics/cmath.h"


ManyBodyOperator OpRadiiAll = {
 name : "RadiiAll",
 rank : 0,
 pi : 0,
 dim : 5,
 size : 5,
 par : NULL,
 me : calcRadii2allod
};


static void ob_radii2all(const SlaterDet* Q,
	  const Gaussian*G1, const Gaussian* G2, 
	  const GaussianAux* X, complex double r2[5])
{
  complex double mr2;
  int A=Q->A, Z=Q->Z, N=Q->N;

  mr2 = (3.0*X->alpha + X->rho2)* X->Q;
 
  r2[0] += 1.0/A*(1.0-1.0/A)* mr2;
  r2[1] += (1.0/(A*A) + (1+G1->xi)/2 *(1.0/Z-2.0/(A*Z)))* mr2;
  r2[2] += (1.0/(A*A) + (1-G1->xi)/2 *(1.0/N-2.0/(A*N)))* mr2;
  r2[3] += 1.0/Z*(1.0-1.0/Z)* (1+G1->xi)/2* mr2;
  r2[4] += 1.0/N*(1.0-1.0/N)* (1-G1->xi)/2* mr2;
}


static void tb_radii2all(const SlaterDet* Q,
	   const Gaussian* G1, const Gaussian* G2, 
	   const Gaussian* G3, const Gaussian* G4, 
	   const GaussianAux* X13, const GaussianAux* X24, 
	   complex double r2[5])
{	
  complex double mr2;
  int A=Q->A, Z=Q->Z, N=Q->N;

  mr2 = cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;

  r2[0] += -2.0/(A*A)* mr2;
  r2[1] +=  2.0*(1.0/(A*A) - 
		 ((1+G1->xi)/2+(1+G2->xi)/2) *1.0/(A*Z))* mr2;
  r2[2] +=  2.0*(1.0/(A*A) - 
		 ((1-G1->xi)/2+(1-G2->xi)/2) *1.0/(A*N))* mr2;
  r2[3] += -2.0/(Z*Z)* (1+G1->xi)/2* (1+G2->xi)/2* mr2;		
  r2[4] += -2.0/(N*N)* (1-G1->xi)/2* (1-G2->xi)/2* mr2;
}


void calcRadii2all(const SlaterDet* Q, const SlaterDetAux* X,
		   RadiiAll* r2)
{
  double r2one[5], r2two[5];
  OneBodyOperator op_ob_radii2 = {dim: 5, opt: 1, par: Q, me: ob_radii2all};
  TwoBodyOperator op_tb_radii2 = {dim: 5, opt: 1, par: Q, me: tb_radii2all};


  calcSlaterDetOBME(Q, X, &op_ob_radii2, r2one);
  calcSlaterDetTBME(Q, X, &op_tb_radii2, r2two);

  r2->r2m = r2one[0] + r2two[0];
  r2->r2p = r2one[1] + r2two[1];
  r2->r2n = r2one[2] + r2two[2];
  r2->r2pp = r2one[3] + r2two[3];
  r2->r2nn = r2one[4] + r2two[4];
}


void calcRadii2allod(void* par,
		     const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, 
		     RadiiAllod* r2)
{
  complex double r2one[5], r2two[5];
  OneBodyOperator op_ob_radii2 = {dim: 5, opt: 1, par: Q, me: ob_radii2all};
  TwoBodyOperator op_tb_radii2 = {dim: 5, opt: 1, par: Q, me: tb_radii2all};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_radii2, r2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_radii2, r2two);

  r2->r2m = r2one[0] + r2two[0];
  r2->r2p = r2one[1] + r2two[1];
  r2->r2n = r2one[2] + r2two[2];
  r2->r2pp = r2one[3] + r2two[3];
  r2->r2nn = r2one[4] + r2two[4];
}


void writeprojectedRadiiAll(FILE* fp,
			    const SlaterDet* Q,
			    const Projection* P,
			    RadiiAllod** radiiall,
			    const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i,n;

  char prefix[8];

  RadiiAllod* r;
  int* idx;
  int ngood;
  complex double *norm, *H;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

        if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+');
        else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+');

        idx = E->index[idxpij(jmax,p,j)];
        norm = E->norm[idxpij(jmax,p,j)];
        H = E->v[idxpij(jmax,p,j)];

	r = radiiall[idxpij(jmax,p,j)];
  
        fprintf(fp, "\n%s  N      = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
        fprintf(fp, "\n%s  H      = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));

        fprintf(fp, "\n%s  Rm     = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(creal(r[idx[i]].r2m/norm[idx[i]])));
        fprintf(fp, "\n%s  Rp     = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(creal(r[idx[i]].r2p/norm[idx[i]])));
        fprintf(fp, "\n%s  Rn     = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(creal(r[idx[i]].r2n/norm[idx[i]])));
        fprintf(fp, "\n%s  Rch    = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(r2charge(creal(r[idx[i]].r2p/norm[idx[i]]), Q->N, Q->Z)));
        fprintf(fp, "\n%s  Rpp    = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(creal(r[idx[i]].r2pp/norm[idx[i]])));
        fprintf(fp, "\n%s  Rnn    = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", sqrt(creal(r[idx[i]].r2nn/norm[idx[i]])));
 
       fprintf(fp, "\n");
      }
    }
}
