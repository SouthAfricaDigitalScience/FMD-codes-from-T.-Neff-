/**

   \file MultiSlaterDet.h

   set of Many-body states described by multiple Slaterdets

   
   (c) 2005 Thomas Neff

*/


#ifndef _MULTISLATERDET_H
#define _MULTISLATERDET_H

#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"
#include "Observables.h"
#include "Symmetry.h"


/// generic description of a set of N many-body states described
/// by a common set of n Slater determinants

#define NMAX 10				/// maximum number of Many-body states per MultiSlaterDet

typedef struct {
  int N;				///< out of N many-body states
  int n;				///< describe n many-body states
  int idx[NMAX];			
} Indices;


struct _MSD_ {
  int A;				///< particle number
  int N;				///< describe N many-body states
  int n;				///< by n SlaterDets
  Symmetry (*symmetry)(const struct _MSD_* MB, int iM);				///< Symmetry 
  complex double (*weight)(const struct _MSD_* MB, int iM, int i);	///< weight of i-th SlaterDet in I-th many-body state
  void (*get)(const struct _MSD_* MB, int i, SlaterDet* Q);		///< copy i-th SlaterDet into Q
  void *internals;			///< for internal use
};

typedef struct _MSD_ MultiSlaterDet;




void extractIndicesfromString(char** str, Indices* Ind);


char* IndicestoStr(const Indices* Ind);


int readMultiSlaterDetfromFile(MultiSlaterDet* Q, Indices* Ind,
			       const char* fname);


void* initprojectedMultiMBME(const Projection* P, const ManyBodyOperator* Op,
			     const MultiSlaterDet* QA, const MultiSlaterDet* QB);

int readprojectedMultiMBMEfromFile(const char* mbfilea, const char* mbfileb,
				   const MultiSlaterDet* QA,
				   const MultiSlaterDet* QB,
				   const Projection* P,
				   const ManyBodyOperator* Op,
				   void** obsme);

void calcprojectedMultiMBME(const Projection* P, const ManyBodyOperator* Op,
			    const MultiSlaterDet* QA, 
			    const MultiSlaterDet* QB,
			    void** obsme);

int writeprojectedMultiMBMEtoFile(const char* mbfilea, const char* mbfileb,
				  const MultiSlaterDet* QA,
				  const MultiSlaterDet* QB,
				  const Projection* P,
				  const ManyBodyOperator* Op,
				  void** obsme);

void calcMultiEigenstatesMulti(const Projection* P,
			       const Interaction* Int,
			       const Observablesod**** obs,
			       const Eigenstates* Ep,
			       const Indices* In,
			       int n,
			       Eigenstates* multiE, Amplitudes* multiA, 
			       double thresh);

void calcexpectprojectedMultiMBME(const Projection* P,
				  const ManyBodyOperator* Op,
				  const void* mbme,
				  const MultiSlaterDet* Q,
				  const Indices* In,
				  int n,
				  const Eigenstates* E,
				  void* expectmbme);

void calctransitionprojectedMultiMBME(const Projection* P,
				      const ManyBodyOperator* Op,
				      const void* mbme,
				      const MultiSlaterDet* Q,
				      const MultiSlaterDet* Qp,
				      const Indices* In,
				      const Indices* Inp,
				      int n, int np,
				      const Eigenstates* E,
				      const Eigenstates* Ep,
				      void* transmbme);

int writeMultiEigenstatestoFile(const char* fname,
				const Projection* P,
				const char* mbfile,
				const Indices* In,
				const Eigenstates* E);

int writeMultiMulticonfigfile(const char* fname,
			      const Projection* P,
			      const char** mbfile,
			      const Indices* In,
			      int n,
			      const Eigenstates* E);

int readMultiMulticonfigfile(const char* fname,
			     char*** mbfilep,
			     Projection* P,
			     MultiSlaterDet** Qp,
			     Indices** Inp,
			     int* nstates,
			     Eigenstates* E);

#endif
