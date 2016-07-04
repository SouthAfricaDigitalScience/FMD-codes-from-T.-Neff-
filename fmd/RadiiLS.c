/**

  \file RadiiLS.c

  calculate spin-orbit contribution to nuclear charge radius


  (c) 2011 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "RadiiLS.h"

#include "misc/physics.h"
#include "numerics/cmath.h"


ManyBodyOperator OpRadiiLS = {
 name : "RadiiLS",
 rank : 0,
 pi : 0,
 dim : 2,
 size : 2,
 par : NULL,
 me : calcRadii2LSod
};


// <R>=0 and <P>=0 assumed implicitly
// calculates l.sigma !

static void ob_radii2ls(void* par,
	  const Gaussian*G1, const Gaussian* G2, 
	  const GaussianAux* X, complex double r2ls[2])
{
  r2ls[0] += (1+G1->xi)/2* cvec3mult(X->rhoxpi, X->sig)* X->T* X->R;
  r2ls[1] += (1-G1->xi)/2* cvec3mult(X->rhoxpi, X->sig)* X->T* X->R;
}


void calcRadii2LS(const SlaterDet* Q, const SlaterDetAux* X,
		  RadiiLS* R2ls)
{
  double r2ls[2];
  OneBodyOperator op_ob_radii2ls = {dim: 2, opt: 0, par: NULL, me: ob_radii2ls};

  calcSlaterDetOBME(Q, X, &op_ob_radii2ls, r2ls);

  R2ls->r2lsp = r2ls[0];
  R2ls->r2lsn = r2ls[1];
}


void calcRadii2LSod(void* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, 
		    RadiiLSod* R2ls)
{
  complex double r2ls[2];
  OneBodyOperator op_ob_radii2ls = {dim: 2, opt: 0, par: NULL, me: ob_radii2ls};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_radii2ls, r2ls);

  R2ls->r2lsp = r2ls[0];
  R2ls->r2lsn = r2ls[1];
}


void writeprojectedRadiiLS(FILE* fp,
			   const SlaterDet* Q,
			   const Projection* P,
			   const RadiiLSod** radiils,
			   const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i,n;

  char prefix[8];

  RadiiLSod* r;
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

	r = radiils[idxpij(jmax,p,j)];
  
        fprintf(fp, "\n%s  N      = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
        fprintf(fp, "\n%s  H      = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));

        fprintf(fp, "\n%s  R2lsp  = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", 0.5*creal(r[idx[i]].r2lsp/norm[idx[i]]));
        fprintf(fp, "\n%s  R2lsn  = ", prefix);
        for (i=0; i<ngood; i++)
          fprintf(fp, "   %8.3f", 0.5*creal(r[idx[i]].r2lsn/norm[idx[i]]));

       fprintf(fp, "\n");
      }
    }
}
