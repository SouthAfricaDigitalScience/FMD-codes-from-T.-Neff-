/**

  \file TimeReversal.c

  calculate time reversal


  (c) 2007 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "TimeReversal.h"

#include "misc/physics.h"


// Time Reversal
ManyBodyOperator OpTimeReversal = {
  name : "TimeRev",
  rank : 0,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcTimeReversalod
};



// these are used as workspace

static SlaterDet Qpp;
static SlaterDetAux Xpp;


void calcTimeReversal(const SlaterDet* Q, const SlaterDetAux* X,
		      complex double* t)
{
  // do we have to initialize the workspace first ?
  if (Qpp.A == 0) {
    initSlaterDet(Q, &Qpp);
    initSlaterDetAux(&Qpp, &Xpp);
  }

  copySlaterDet(Q, &Qpp);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  double norm = Xpp.ovlap;

  timerevertSlaterDet(&Qpp);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  *t = Xpp.ovlap/norm;
}


void calcTimeReversalod(void* dummy,
			const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X, 
			complex double* t)
{
  // do we have to initialize the workspace first ?
  if (Qpp.A == 0) {
    initSlaterDet(Qp, &Qpp);
    initSlaterDetAux(&Qpp, &Xpp);
  }

  copySlaterDet(Qp, &Qpp);
  timerevertSlaterDet(&Qpp);

  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  *t = Xpp.ovlap;
}



void showprojectedTimeReversal(FILE* fp,
			       const Projection* P,
			       complex double** tr,
			       const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,ipj,i;

  char prefix[20];

  complex double* norm;
  complex double* H;
  complex double* TR;
  int* idx;
  int ngood;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      ipj = idxpij(jmax,p,j);
      ngood = E->ngood[ipj];

      if (ngood) {

        if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
        else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

        idx = E->index[idxpij(jmax,p,j)];
        norm = E->norm[idxpij(jmax,p,j)];
        H = E->v[idxpij(jmax,p,j)];
 
        TR = tr[idxpij(jmax,p,j)];

        fprintf(fp, "\n%s N     = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
        fprintf(fp, "\n%s H     = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
        fprintf(fp, "\n%s TR    = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   (%8.3f, %8.3f)",
		  creal(TR[idx[i]]/norm[idx[i]]),
		  cimag(TR[idx[i]]/norm[idx[i]]));

        fprintf(fp, "\n");      
      }
    }
}
