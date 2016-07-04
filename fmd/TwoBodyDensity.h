/**

  \file TwoBodyDensity.h

  calculate matrix elements of twob-body density operators


  (c) 20010 Thomas Neff

*/


#ifndef _TWOBODYDENSITY_H
#define _TWOBODYDENSITY_H


#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"

extern ManyBodyOperator OpPairs;


typedef struct {
  double rmax;
  int npoints;
} TBDensRPara;

extern ManyBodyOperator OpTwoBodyDensityR;

typedef struct {
  double qmax;
  int npoints;
} TBDensQPara;

extern ManyBodyOperator OpTwoBodyDensityQ;


typedef struct {
  double rmax;
  int npoints;
  int lmax;
  int lambdamax;
} TBDensRLPara;

extern ManyBodyOperator OpTwoBodyDensityRL;


typedef struct {
  double qmax;
  int npoints;
  int lmax;
  int lambdamax;
} TBDensQLPara;

extern ManyBodyOperator OpTwoBodyDensityQL;



void calcPairs(const SlaterDet* Q, const SlaterDetAux* X, double pairs[4]);

void calcPairsod(void* Par,
		 const SlaterDet* Q, const SlaterDet* Qp,
		 const SlaterDetAux* X, complex double pairs[4]);

void writeprojectedPairs(FILE* fp, const Projection* P, 
			 const complex double (**pairs)[4],
			 const Eigenstates* E);

void initOpTwoBodyDensityR(TBDensRPara* par);

void calcTBDensRod(TBDensRPara* par,
                   const SlaterDet* Q, const SlaterDet* Qp,
                   const SlaterDetAux* X, complex double* densr);

void writeTBDensR(FILE* fp,
		  const Projection* P,
		  const TBDensRPara* p,
		  int j, int pi, int a,
		  void* tbdensme,
		  const Eigenstates* E);


void initOpTwoBodyDensityQ(TBDensQPara* par);

void calcTBDensQod(TBDensQPara* par,
                   const SlaterDet* Q, const SlaterDet* Qp,
                   const SlaterDetAux* X, complex double* densq);

void writeTBDensQ(FILE* fp,
		  const Projection* P,
		  const TBDensQPara* p,
		  int j, int pi, int a,
		  void* tbdensme,
		  const Eigenstates* E);


void initOpTwoBodyDensityRL(TBDensRLPara* par);

void calcTBDensRLod(TBDensRLPara* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, complex double* densr);

void writeTBDensRL(FILE* fp,
		   const Projection* P,
		   const TBDensRLPara* p,
		   int j, int pi, int a,
		   void* tbdensme,
		   const Eigenstates* E);


void initOpTwoBodyDensityQL(TBDensQLPara* par);

void calcTBDensQLod(TBDensQLPara* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, complex double* densq);

void writeTBDensQL(FILE* fp,
		   const Projection* P,
		   const TBDensQLPara* p,
		   int j, int pi, int a,
		   void* tbdensme,
		   const Eigenstates* E);

#endif
