/**

  \file ProjectedObservables.h

  read or calculate/write projected Observables matrix elements


  (c) 2003 Thomas Neff

*/


#ifndef _PROJECTEDOBSERVABLES_H
#define _PROJECTEDOBSERVABLES_H

#include "SlaterDet.h"
#include "Observables.h"
#include "Projection.h"
#include "Symmetry.h"


extern ManyBodyOperator OpObservables;


void initOpObservables(const Interaction* Int);


void scaleprojectedObservablesMBME(const Projection* P,
                                   const Interaction* Int,
                                   Observablesod** obsme);


void showprojectedObservables(FILE* fp,
			      const Projection* P,
			      const Interaction* Int,
			      const SlaterDet* Q,
			      const Observablesod* obs,
			      const Eigenstates* E,
			      const Amplitudes* A,
			      const char* pre);

void showSpectrum(FILE* fp,
                  const Projection* P,
                  const Observablesod** obs,
                  const Eigenstates* E);

#endif
