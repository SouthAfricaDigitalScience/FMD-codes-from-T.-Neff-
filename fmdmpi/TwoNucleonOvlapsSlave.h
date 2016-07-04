/**

  \file TwoNucleonOvlapsSlave.h

  MPI slave for calculating Two-Nucleon Ovlap MEs 


  (c) 2012 Thomas Neff

*/


#ifndef _TWONUCLEONOVLAPSSLAVE_H
#define _TWONUCLEONOVLAPSSLAVE_H

#include "fmd/TwoNucleonOvlaps.h"


void BroadcastTwoNucleonOvlapsPara(TwoNucleonOvlapsPara* TNOpar);

void TwoNucleonOvlapsSlave(void);


#endif
