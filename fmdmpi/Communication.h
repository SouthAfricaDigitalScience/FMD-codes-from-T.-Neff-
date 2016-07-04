/**

  \file Communication.h

  Communicate using MPI


  (c) 2003 Thomas Neff

*/


#ifndef _COMMUNICATION_H
#define _COMMUNICATION_H


#define TASKFIN -1
#define TASKSTART 1
#define TASKSTARTCM 2
#define TASKPOTENTIAL 10
#define TASKGRADPOTENTIAL 11
#define TASKPOTENTIALOD 12
#define TASKGRADPOTENTIALOD 13
#define TASKHAMILTONIANOD 14
#define TASKGRADHAMILTONIANOD 15

#define TASKPROJECTOVLAPOD 100
#define TASKPROJECTOBSERVABLESOD 101
#define TASKPROJECTMONOPOLEFORMFACTOROD 110
#define TASKPROJECTEDDIAGONALDENSITYMATRIXHOOD 120
#define TASKPROJECTONENUCLEONOVLAPSOD 130
#define TASKPROJECTTWONUCLEONOVLAPSTOD 140
#define TASKPROJECTTWONUCLEONOVLAPSYOD 141


#define TAGROW 0
#define TAGROWCOL 1
#define TAGANGLES3 2
#define TAGPROJECT3 3
#define TAGPROJECT6 4
#define TAGPOTENTIAL 11
#define TAGGRADPOTENTIALVAL 12
#define TAGGRADPOTENTIALGRAD 13
#define TAGPOTENTIALOD 21
#define TAGGRADPOTENTIALODVAL 22
#define TAGGRADPOTENTIALODGRAD 23
#define TAGOVLAPOD 31
#define TAGGRADOVLAPODVAL 32
#define TAGGRADOVLAPODGRAD 33
#define TAGHAMILTONIANOD 41
#define TAGGRADHAMILTONIANODVAL 42
#define TAGGRADHAMILTONIANODGRAD 43
#define TAGANGULARMOMENTUMOD 51 

#define TAGMEOD 1001


#include "../fmd/Interaction.h"
#include "../fmd/SlaterDet.h"
#include "../fmd/gradSlaterDet.h"


extern int mpirank;
extern int mpisize;


/// use MPI Broadcast to distribute Interaction parameters
/// slaves allocate space for Interaction
/// labels are not broadcastet
void BroadcastInteraction(Interaction* Int);

/// generic broadcast for Parameters
void BroadcastParameters(void* par, int size);

/// broadcast size of nucleus
void BroadcastA(int* A);

void BroadcastTask(int* task);

// slaves overwride sldet
void BroadcastSlaterDet(SlaterDet* Q);

void BroadcastSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X);

void BroadcastgradSlaterDetAux(const SlaterDet* Q, gradSlaterDetAux* dX);
  
#endif

