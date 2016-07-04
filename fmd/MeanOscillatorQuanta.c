/**

   \file MeanOscillatorQuanta.c

   calculate mean number of oscillator quanta in many-body state


   (c) 2007 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "MeanOscillatorQuanta.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


// r2 and p2 operators
ManyBodyOperator OpMeanOsciQuanta = {
 name : "MeanOsciQuanta",
 rank : 0,
 pi : 0,
 dim : 6,
 size : 6,
 par : NULL,
 me : calcMeanOsciQuantaod
};


static void ob_moscquanta(void* par,
			  const Gaussian *G1, const Gaussian* G2,
			  const GaussianAux *X, MeanOsciQuanta* q)
{
  complex double mr2, mp2;

  mr2 = (3.0*X->alpha + X->rho2)* X->Q;
  mp2 = (3*X->lambda + X->pi2) * X->Q;

  q->r2[0] = mr2;
  q->r2[1] = (1+G1->xi)/2* mr2;
  q->r2[2] = (1-G1->xi)/2* mr2;

  q->p2[0] = mp2;
  q->p2[1] = (1+G1->xi)/2* mp2;
  q->p2[2] = (1-G1->xi)/2* mp2;
}


void calcMeanOsciQuantaod(void* Par,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X,
			  MeanOsciQuanta* q)
{
  OneBodyOperator op_ob_moscquanta = {dim: 6, opt: 1, par: NULL, me: ob_moscquanta};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_moscquanta, q);
}
			      

static inline double dsqr(double x)
{
  return x*x;
}


void writeprojectedMeanOscillatorQuanta(FILE* fp, double omega,
					const Projection* P,
					const MeanOsciQuanta** Q,
					const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i;

  char prefix[8];

  int* idx;
  int ngood;
  complex double *norm, *H;
  MeanOsciQuanta *q;
  

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

	if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
	else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
	norm = E->norm[idxpij(jmax,p,j)];
	H = E->v[idxpij(jmax,p,j)];

	q = Q[idxpij(jmax,p,j)];

	fprintf(fp, "\n%s N     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
	fprintf(fp, "\n%s H     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
	fprintf(fp, "\n%s Q     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", 
		  0.5*(
		       creal(q[idx[i]].p2[0]/norm[idx[i]])/(mnucleon*omega) +
		       mnucleon*omega*creal(q[idx[i]].r2[0]/norm[idx[i]])
		       ));
	fprintf(fp, "\n%s Qp    = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", 
		  0.5*(
		       creal(q[idx[i]].p2[1]/norm[idx[i]])/(mnucleon*omega) +
		       mnucleon*omega*creal(q[idx[i]].r2[1]/norm[idx[i]])
		       ));
	fprintf(fp, "\n%s Qn    = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", 
		  0.5*(
		       creal(q[idx[i]].p2[2]/norm[idx[i]])/(mnucleon*omega) +
		       mnucleon*omega*creal(q[idx[i]].r2[2]/norm[idx[i]])
		       ));

	fprintf(fp, "\n");	  
      }
    }
}
