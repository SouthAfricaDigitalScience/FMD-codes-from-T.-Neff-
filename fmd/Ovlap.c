/**

  \file Overlap.c

  calculate Overlap for Slater determinants


  (c) 2003 Thomas Neff

*/

#include <math.h>

#include "SlaterDet.h"
#include "Projection.h"
#include "Ovlap.h"
#include "misc/physics.h"


// Ovlap
ManyBodyOperator OpOvlap = {
  name : "Ovlap",
  rank : 0,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcOvlapod
};


// extract Ovlap from SlaterDetAux
void calcOvlapod(void* dummy,
		 const SlaterDet* Q, const SlaterDet* Qp,
		 const SlaterDetAux* X,
		 complex double* ovl)
{
  *ovl = X->ovlap;
}


void writeprojectedtransitionOvlap(FILE* fp,
                                   const Projection* P,
                                   const complex double**** ovlme,
                                   const Eigenstates* Efin, 
                                   const Eigenstates* Eini)
{
  int odd = P->odd;
  int jmax = P->jmax;

  int p, j;
  int idx;
  int *idxi, *idxf;
  int i, f;
  complex double O;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      idx = idxpij(jmax,p,j);        

      if (Eini->ngood[idx] == 0 || Efin->ngood[idx] == 0)
	    
        continue;

      idxi = Eini->index[idx];
      idxf = Efin->index[idx];

      if (P->odd) fprintf(fp, "\n\n Ovlap [%d/2%c]\n\n",
                          j, p ? '-' : '+'); 
      else	  fprintf(fp, "\n\n Ovlap [%d%c]\n\n",
                          j/2, p ? '-' : '+'); 

      fprintf(fp, "%14c", ' ');
      for (i=0; i<Eini->ngood[idx]; i++) 
        fprintf(fp, "%10.3f MeV", hbc*creal(Eini->v[idx][idxi[i]]));

      for (f=0; f<Efin->ngood[idx]; f++) {
        fprintf(fp, "\n%10.3f MeV", hbc*creal(Efin->v[idx][idxf[f]]));
        for (i=0; i<Eini->ngood[idx]; i++) {
          O = ovlme[idx][idx][idxf[f]][idxi[i]]/
            sqrt(cabs(Efin->norm[idx][idxf[f]]*Eini->norm[idx][idxi[i]]));
          fprintf(fp, "%10.5f    ", cabs(O));
	}
	fprintf(fp, "\n");	
      }
    }
}
