/**

   \file SymmetricSlaterDet.c

   implement SymmetricSlaterDet as MultiSlaterDet

   only parity for the start
   
   (c) 2005

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"
#include "SymmetricSlaterDet.h"

#include "misc/utils.h"


typedef struct {
  SlaterDet* Q;
  Symmetry S;
} SymmetricSlaterDetInternals;


/// some fields have to be initialized later

MultiSlaterDet SymmetricSlaterDet = {
  A : 0,
  N : 1,
  n : 1,
  symmetry : SymmetricSlaterDetSymmetry,
  weight : SymmetricSlaterDetWeight,
  get : SymmetricSlaterDetGet,
  internals : NULL			// pointer to SymmetricSlaterDetInternals
};


Symmetry SymmetricSlaterDetSymmetry(const MultiSlaterDet* MB, int iM)
{
  SymmetricSlaterDetInternals* internals = MB->internals;

  return internals->S;
}


complex double SymmetricSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i)
{
  return 1.0;
}


void SymmetricSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q)
{
  SymmetricSlaterDetInternals* internals = MB->internals;

  copySlaterDet(internals->Q, Q);
}


int SymmetricSlaterDetRead(FILE* fp, MultiSlaterDet* MB)
{
  *MB = SymmetricSlaterDet;

  SymmetricSlaterDetInternals* internals = 
    malloc(sizeof(SymmetricSlaterDetInternals));

  internals->Q = malloc(sizeof(SlaterDet));

  char buf[1024];
  fgets(buf, 1024, fp);
  char sym[255]; char* s=sym;
  sscanf(buf, "<SymmetricSlaterDet Sym=%s>", sym);
  extractSymmetryfromString(&s, &internals->S);

  int ret;
  ret = readSlaterDet(fp, internals->Q);

  MB->A = internals->Q->A;
  MB->internals = internals;

  return ret;
}
