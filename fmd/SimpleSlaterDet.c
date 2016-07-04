/**

   \file SimpleSlaterDet.c

   implement SlaterDet as MultiSlaterDet

   this is only glue
   
   (c) 2005

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"
#include "MultiSlaterDet.h"
#include "SimpleSlaterDet.h"


typedef struct {
  SlaterDet* Q;
} SimpleSlaterDetInternals;


/// some fields have to be initialized later

MultiSlaterDet SimpleSlaterDet = {
  A : 0,
  N : 1,
  n : 1,
  symmetry : SimpleSlaterDetSymmetry,
  weight : SimpleSlaterDetWeight,
  get : SimpleSlaterDetGet,
  internals : NULL			// pointer to SimpleSlaterDetInternals
};


Symmetry SimpleSlaterDetSymmetry(const MultiSlaterDet* MB, int iM)
{
  return 0;
}


complex double SimpleSlaterDetWeight(const MultiSlaterDet* MB, int iM, int i)
{
  return 1.0;
}


void SimpleSlaterDetGet(const MultiSlaterDet* MB, int i, SlaterDet* Q)
{
  SimpleSlaterDetInternals* internals = MB->internals;

  copySlaterDet(internals->Q, Q);
}



int SimpleSlaterDetRead(FILE* fp, MultiSlaterDet* MB)
{
  *MB = SimpleSlaterDet;

  SimpleSlaterDetInternals* internals = 
    malloc(sizeof(SimpleSlaterDetInternals));

  internals->Q = malloc(sizeof(SlaterDet));

  MB->internals = internals;

  int ret;
  ret = readSlaterDet(fp, internals->Q);

  MB->A = internals->Q->A;

  return ret;
}
