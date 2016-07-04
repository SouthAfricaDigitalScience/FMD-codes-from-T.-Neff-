/**

  \file GamovTeller.c 

  GamovTeller transitions

  calculates Gamov Teller transition matrix elements
  in spherical basis


  (c) 2004,2007 Thomas Neff

*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "GamovTeller.h"

#include "numerics/cmath.h"
#include "misc/utils.h"
#include "misc/physics.h"


// GTplus sigma tau+
ManyBodyOperator OpGTplus = {
  name : "GTplus",
  rank : 2,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcGTplusod
};

// GTminus sigma tau-
ManyBodyOperator OpGTminus = {
  name : "GTminus",
  rank : 2,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcGTminusod
};


static void ob_gtplus(void* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, complex double gtp[3])
{
  int i;

  if (G1->xi == 1 && G2->xi == -1)
    for (i=0; i<3; i++)
      gtp[i] += X->sig[i]*X->R;
}


static void ob_gtminus(void* par,
		       const Gaussian* G1, const Gaussian* G2, 
		       const GaussianAux* X, complex double gtm[3])
{
  int i;

  if (G1->xi == -1 && G2->xi == 1)
    for (i=0; i<3; i++)
      gtm[i] += X->sig[i]*X->R;
}



void calcGTplusod(void* Par,
		  const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double GTplus[3])
{
  complex double gtp[3];

  OneBodyOperator op_ob_gtp = {dim: 3, opt: 0, par: NULL, me: ob_gtplus};

  // SlaterDetAux X calculated by calcSlaterDetAuxodSingular needed here
  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_gtp, gtp);

  // spherical components -1, 0, +1
  GTplus[0] = 	sqrt(0.5)*(gtp[0] - I*gtp[1]);
  GTplus[1] =   gtp[2];
  GTplus[2] = - sqrt(0.5)*(gtp[0] + I*gtp[1]);
}


void calcGTminusod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double GTminus[3])
{
  complex double gtm[3];

  OneBodyOperator op_ob_gtm = {dim: 3, opt: 0, par: NULL, me: ob_gtminus};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_gtm, gtm);

  // spherical components -1, 0, +1
  GTminus[0] = 	 sqrt(0.5)*(gtm[0] - I*gtm[1]);
  GTminus[1] =   gtm[2];
  GTminus[2] = - sqrt(0.5)*(gtm[0] + I*gtm[1]);
}

					 

#define SQR(x) ((x)*(x))


void writeprojectedtransitions(FILE* fp,
			       const Projection* P,
			       int rank, int pi,
			       const char* label, const char* unit,
			       const complex double**** transme,
			       const Eigenstates* Efin, 
			       const Eigenstates* Eini)
{
  int odd = P->odd;
  int jmax = P->jmax;

  int pini, jini, pfin, jfin;
  int idxini, idxfin;
  int *idxi, *idxf;
  int i, f;
  double B;

  for (pini=0; pini<=1; pini++)
    for (jini=odd; jini<jmax; jini=jini+2) {
      pfin = (pini+pi)%2;
      for (jfin=abs(jini-rank); jfin<=min(jini+rank,jmax-1); jfin=jfin+2) {

	idxini = idxpij(jmax,pini,jini);
	idxfin = idxpij(jmax,pfin,jfin);

	if (Eini->ngood[idxini] == 0 || Efin->ngood[idxfin] == 0)
	    
	  continue;

	idxi = Eini->index[idxini];
	idxf = Efin->index[idxfin];

	if (P->odd) fprintf(fp, "\n\n %s [%d/2%c] <- [%d/2%c] (%s)\n\n",
			    label, jfin, pfin ? '-' : '+', 
			    jini, pini ? '-' : '+', unit);
	else	    fprintf(fp, "\n\n %s [%d%c] <- [%d%c] (%s)\n\n",
			    label, jfin/2, pfin ? '-' : '+', 
			    jini/2, pini ? '-' : '+', unit);

	fprintf(fp, "%14c", ' ');
	for (i=0; i<Eini->ngood[idxini]; i++) 
	  fprintf(fp, "%10.3f MeV", hbc*creal(Eini->v[idxini][idxi[i]]));

	for (f=0; f<Efin->ngood[idxfin]; f++) {
	  fprintf(fp, "\n%10.3f MeV", hbc*creal(Efin->v[idxfin][idxf[f]]));
	  for (i=0; i<Eini->ngood[idxini]; i++) {
            B = SQR(gagv)*(jfin+1.0)/(jini+1.0)*
              SQR(cabs(transme[idxfin][idxini][idxf[f]][idxi[i]]))/
              cabs(Efin->norm[idxfin][idxf[f]]*Eini->norm[idxini][idxi[i]]);
	    if (creal(Efin->v[idxfin][idxf[f]]) < creal(Eini->v[idxini][idxi[i]]))
	      fprintf(fp, "%10.5f    ", B);
	    else
	      fprintf(fp, "%10.5f *  ", B);
	  }
	}
	fprintf(fp, "\n");	
      }
    }
}


void writeprojectedtransitionGTplus(FILE* fp,
				    const Projection* P,
				    const complex double**** gtp,
				    const Eigenstates* Efin, 
				    const Eigenstates* Eini)
{
  fprintf(fp, "\n ########### Gamov Teller plus transitions\n");
  writeprojectedtransitions(fp, P, 2, 0, "B(GT+)", "", gtp, Efin, Eini);
}

void writeprojectedtransitionGTminus(FILE* fp,
                                     const Projection* P,
                                     const complex double**** gtm,
                                     const Eigenstates* Efin, 
                                     const Eigenstates* Eini)
{
  fprintf(fp, "\n ########### Gamov Teller minus transitions\n");
  writeprojectedtransitions(fp, P, 2, 0, "B(GT-)", "", gtm, Efin, Eini);
}


