/**

  \file Radii.c

  calculate radii without correcting for center of mass


  (c) 2003 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "SDRadii.h"

#include "misc/physics.h"
#include "numerics/cmath.h"


ManyBodyOperator OpSDRadii = {
 name : "SDRadii",
 rank : 0,
 pi : 0,
 dim : 3,
 size : 3,
 par : NULL,
 me : calcSDRadii2od
};


// do not correct for center of mass
static void ob_radii2sd(const SlaterDet* Q,
	  const Gaussian*G1, const Gaussian* G2, 
	  const GaussianAux* X, complex double r2[3])
{
  complex double mr2;
  int A=Q->A, Z=Q->Z, N=Q->N;

  mr2 = (3.0*X->alpha + X->rho2)* X->Q;
 
  r2[0] += 1.0/A* mr2;
  r2[1] += 1.0/Z* (1+G1->xi)/2*  mr2;
  r2[2] += 1.0/N* (1-G1->xi)/2* mr2;
}


void calcSDRadii2(const SlaterDet* Q, const SlaterDetAux* X, 
		  SDRadii* r2)
{
  double r2sd[3];
  OneBodyOperator op_ob_radii2sd = {dim: 3, opt: 1, par: Q, me: ob_radii2sd};


  calcSlaterDetOBME(Q, X, &op_ob_radii2sd, r2sd);

  r2->r2m = r2sd[0];
  r2->r2p = r2sd[1];
  r2->r2n = r2sd[2];
}


void calcSDRadii2od(void* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X,
		    SDRadiiod* r2)
{
  complex double r2sd[3];
  OneBodyOperator op_ob_radii2sd = {dim: 3, opt: 1, par: Q, me: ob_radii2sd};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_radii2sd, r2sd);

  r2->r2m = r2sd[0];
  r2->r2p = r2sd[1];
  r2->r2n = r2sd[2];
}


void writeprojectedSDRadii(FILE* fp,
			   const SlaterDet* Q,
			   const Projection* P,
			   SDRadiiod** sdradii,
			    const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i,n;

  char prefix[8];

  SDRadiiod* r;
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

	r = sdradii[idxpij(jmax,p,j)];
  
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
 
       fprintf(fp, "\n");
      }
    }
}
