/**

  \file ProjectedObservables.c

  read or calculate/write projected Observables matrix elements


  (c) 2003 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "SlaterDet.h"
#include "Observables.h"
#include "Projection.h"
#include "Symmetry.h"
#include "Interaction.h"

#include "misc/utils.h"
#include "misc/physics.h"


// has to be initialized with initOpObservables

ManyBodyOperator OpObservables = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcObservablesod
};


void initOpObservables(const Interaction* Int)
{
  char* obsintname = malloc(13+strlen(Int->name));
  sprintf(obsintname, "Observables-%s", Int->name);
  OpObservables.name = obsintname;
  OpObservables.dim = sizeof(Observablesod)/sizeof(double complex) - (MAXINTERACTIONS+1) + Int->n ;
  OpObservables.size = sizeof(Observablesod)/sizeof(double complex);
  OpObservables.par = Int;
}


void scaleprojectedObservablesMBME(const Projection* P,
                                   const Interaction* Int,
                                   Observablesod** obsme)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,m,k;
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2)
      for (m=-j; m<=j; m=m+2)
        for (k=-j; k<=j; k=k+2)
          scaleObservablesod(Int, &obsme[idxpij(jmax,p,j)][idxjmk(j,m,k)]);
}


void showprojectedObservables(FILE* fp,
			      const Projection* P,
			      const Interaction* Int,
			      const SlaterDet* Q,
			      const Observablesod** obs,
			      const Eigenstates* E,
			      const Amplitudes* A,
			      const char* pre)
{
  int odd=P->odd;
  int jmax=P->jmax;
  int n=E->n;

  int p,j,ipj,i,a,ai;

  char prefix[20];

  Observablesod* o;
  int* idx;
  int ngood, ngooda;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      ipj = idxpij(jmax,p,j);
      ngood = E->ngood[ipj];

      if (ngood) {

	if(odd) sprintf(prefix, "%s[%d/2%c]", pre, j, p ? '-' : '+'); 
	else    sprintf(prefix, "%s[%d%c]", pre, j/2, p ? '-' : '+'); 

	o = obs[ipj]; 
	idx = E->index[ipj];

	fprintf(fp, "\n%s  N         = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.6f", creal(o[idx[i]].n));
	fprintf(fp, "\n%s  H         = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(o[idx[i]].h/o[idx[i]].n));
	fprintf(fp, "\n%s  Tcm       = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(o[idx[i]].tcm/o[idx[i]].n));
	fprintf(fp, "\n%s  T         = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(o[idx[i]].t/o[idx[i]].n));
	fprintf(fp, "\n%s  sumV      = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(o[idx[i]].v[0]/o[idx[i]].n));
	for (a=1; a<Int->n; a++) {
	  fprintf(fp, "\n%s  %-10s= ", prefix, Int->label[a]); 
	  for (i=0; i<ngood; i++)
	    fprintf(fp, "   %8.3f", hbc*creal(o[idx[i]].v[a]/o[idx[i]].n));
	}
	fprintf(fp, "\n%s  L2        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(o[idx[i]].l2/o[idx[i]].n));
	fprintf(fp, "\n%s  S2        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(o[idx[i]].s2/o[idx[i]].n));
	fprintf(fp, "\n%s  J2        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(o[idx[i]].j2/o[idx[i]].n));
	fprintf(fp, "\n%s  Pi        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(o[idx[i]].pi/o[idx[i]].n));
	fprintf(fp, "\n%s  T2        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(o[idx[i]].t2/o[idx[i]].n));
	fprintf(fp, "\n%s  Rm        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", sqrt(creal(o[idx[i]].r2m/o[idx[i]].n)));
	fprintf(fp, "\n%s  Rp        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", sqrt(creal(o[idx[i]].r2p/o[idx[i]].n)));
	fprintf(fp, "\n%s  Rn        = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", sqrt(creal(o[idx[i]].r2n/o[idx[i]].n)));
	fprintf(fp, "\n%s  Rch       = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", sqrt(r2charge(creal(o[idx[i]].r2p/o[idx[i]].n), Q->N, Q->Z)));

	// Amplitudes ?
	if (A) {
	  for (a=0; a<n; a++)
	    if (ngooda = A->ngood[ipj][a]) {
	      for (ai=0; ai<ngooda; ai++) {
		fprintf(fp, "\n%s  a[%2d][%d]  = ", prefix, a, ai);
		for (i=0; i<ngood; i++)
		  fprintf(fp, "   %8.3f", 
			  cabs(A->amp[ipj][ai+a*(j+1)+idx[i]*n*(j+1)]));
	      }
	    }
	}

	fprintf(fp, "\n");	  
      }
    }
}


typedef struct {
  double en;
  char spin[6];
  int i;
  double j2;
  double l2;
  double s2;
  double t2;
  double tcm;
} Level;

static int cmplevels(Level* a, Level *b)
{
  return (a->en >= b->en ? (a->en > b->en ? 1 : 0) : -1);
}

void showSpectrum(FILE* fp,
                  const Projection* P,
                  const Observablesod** obs,
                  const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,ipj,i;

  char spin[6];

  Observablesod* o;
  int* idx;

  // how many levels ?
  int nlevels = 0;
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2)
      nlevels += E->ngood[idxpij(jmax,p,j)];

  Level level[nlevels];

  // collect level data
  int ilevel = -1;
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      ipj = idxpij(jmax,p,j);
      if(odd) sprintf(spin, "%d/2%c", j, p ? '-' : '+'); 
      else    sprintf(spin, "%d%c", j/2, p ? '-' : '+'); 

      o = obs[ipj]; 
      idx = E->index[ipj];

      for (i=0; i<E->ngood[ipj]; i++) {
        ilevel++;

        level[ilevel].en = hbc*creal(o[idx[i]].h/o[idx[i]].n);
        strncpy(level[ilevel].spin, spin, 6); 
        level[ilevel].i = i;
        level[ilevel].j2 = creal(o[idx[i]].j2/o[idx[i]].n);
        level[ilevel].l2 = creal(o[idx[i]].l2/o[idx[i]].n);
        level[ilevel].s2 = creal(o[idx[i]].s2/o[idx[i]].n);
        level[ilevel].t2 = creal(o[idx[i]].t2/o[idx[i]].n);
        level[ilevel].tcm = hbc*creal(o[idx[i]].tcm/o[idx[i]].n);
      }
    }

  // sort levels according to energy
  qsort(level, nlevels, sizeof(Level), cmplevels);

  double e0 = level[0].en;

  // write level data

  fprintf(fp, "\n# E0 = %8.3f MeV\n", e0);
  fprintf(fp, "# Energy\t  Spin   #\t L2\t S2\t T2\t\t J2\t Tcm\n\n");

  for (i=0; i<nlevels; i++) {
    fprintf(fp, "%8.3f\t%6s  %2d\t%6.3f\t%6.3f\t%6.3f\t\t%6.3f\t%7.3f\n",
            level[i].en - e0, level[i].spin, level[i].i, level[i].l2, level[i].s2, level[i].t2,
            level[i].j2, level[i].tcm);
  } 

}
  
