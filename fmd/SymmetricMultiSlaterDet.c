/**

   \file SymmetricMultiSlaterDet.c

   implement multiconfig state as MultiSlaterDet


   (c) 2006 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"
#include "SymmetricMultiSlaterDet.h"

#include "misc/utils.h"


typedef struct {
  SlaterDet* Q;
  complex double* w;
  Symmetry S;
} SymmetricMultiSlaterDetInternals;


/// some fields have to be initialized later

MultiSlaterDet SymmetricMultiSlaterDet = {
  A : 0,
  N : 1,
  n : 0,
  symmetry : SymmetricMultiSlaterDetSymmetry,
  weight : SymmetricMultiSlaterDetWeight,
  get : SymmetricMultiSlaterDetGet,
  internals : NULL			// pointer to SymmetricMultiSlaterDetInternals
};


Symmetry SymmetricMultiSlaterDetSymmetry(const MultiSlaterDet* MB, int iM)
{
  SymmetricMultiSlaterDetInternals* internals = MB->internals;

  return internals->S;
}


complex double SymmetricMultiSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i)
{
  SymmetricMultiSlaterDetInternals* internals = MB->internals;

  return internals->w[i];
}


void SymmetricMultiSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q)
{
  SymmetricMultiSlaterDetInternals* internals = MB->internals;

  copySlaterDet(&internals->Q[i], Q);
}


int SymmetricMultiSlaterDetRead(FILE* fp, MultiSlaterDet* MB)
{
  *MB = SymmetricMultiSlaterDet;

  SymmetricMultiSlaterDetInternals* internals = 
    malloc(sizeof(SymmetricMultiSlaterDetInternals));

  char buf[1024];
  fgets(buf, 1024, fp);

  int n;
  char sym[255]; char* s=sym;
  sscanf(buf, "<SymmetricMultiSlaterDet n=%d Sym=%s>", &n, sym);
  extractSymmetryfromString(&s, &internals->S);

  MB->n = n;
  internals->Q = malloc(n*sizeof(SlaterDet));
  internals->w = malloc(n*sizeof(complex double));

  do
    fgets(buf, 1024, fp);
  while (strncmp(buf, "<Weights", 8));
  int i;
  double wre, wim;
  char* c = strtok(buf+8, ", ()");
  for (i=0; i<n; i++) {
    wre=atof(c); c=strtok(NULL, " ,()");
    wim=atof(c); c=strtok(NULL, " ,()");
    internals->w[i] = wre+I*wim;
  }

  int ret=0;
  for (i=0; i<n; i++)
    ret |= readSlaterDet(fp, &internals->Q[i]);

  MB->A = internals->Q[0].A;
  MB->internals = internals;

  return ret;
}
