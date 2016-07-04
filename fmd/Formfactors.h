/**

   \file Formfactors.h

   calculate formfactors


   (c) 2004-2007 Thomas Neff

*/


#ifndef _FORMFACTORS_H
#define _FORMFACTORS_H

#include <stdio.h>

#include "SlaterDet.h"
#include "Projection.h"


#define NMULTIPOLES 4


typedef struct {
  double qmax;
  int npoints;
  int nalpha;
  int ncosb;
  int recoil;
} FormfactorPara;


extern ManyBodyOperator OpMultipoleFormfactor[];

extern ManyBodyOperators OpMultipoleFormfactors;



void initOpFormfactors(FormfactorPara* par);


void calcMultipoleFormfactorsod(FormfactorPara* par,
				const SlaterDet* Q, const SlaterDet* Qp,
				const SlaterDetAux* X,
				complex double* ffactor);


void writeChargeFormfactors(FILE* fp,
			    const Projection* P,
			    const FormfactorPara* par,
			    int l,
			    int j, int pi, int a,
			    void* ffactorme,
			    const Eigenstates* E);

void writeTransitionChargeFormfactors(FILE* fp,
				      const Projection* P,
				      const FormfactorPara* p,
				      int l,
				      int jfin, int pfin, int afin,
				      int jini, int pini, int aini,
				      void* ffactorme,
				      const Eigenstates* E);

void writePointFormfactors(FILE* fp,
			   const Projection* P,
			   const FormfactorPara* par,
			   int l,
			   int j, int pi, int a,
			   void* ffactorme,
			   const Eigenstates* E);

void writeTransitionPointFormfactors(FILE* fp,
				      const Projection* P,
				      const FormfactorPara* p,
				      int l,
				      int jfin, int pfin, int afin,
				      int jini, int pini, int aini,
				      void* ffactorme,
				      const Eigenstates* E);

void writeChargeDensities(FILE* fp,
			  const Projection* P,
			  const FormfactorPara* p,
			  double rmax, int npointsr,
			  int A, int l,
			  int j, int pi, int a,
			  void* ffactorme,
			  const Eigenstates* E);

void writeChargeFourierBessel(FILE* fp,
			      const Projection* P,
			      const FormfactorPara* p,
			      double rmax, int npointsr,
			      int A, int l,
			      int j, int pi, int a,
			      void* ffactorme,
			      const Eigenstates* E);

void writePointDensities(FILE* fp,
			 const Projection* P,
			 const FormfactorPara* p,
			 double rmax, int npointsr,
			 int A, int l,
			 int j, int pi, int a,
			 void* ffactorme,
			 const Eigenstates* E);

void writeMatterDensities(FILE* fp,
                          const Projection* P,
                          const FormfactorPara* p,
                          double rmax, int npointsr,
                          int A, int l,
                          int j, int pi, int a,
                          void* ffactorme,
                          const Eigenstates* E);

void writeTransitionPointDensities(FILE* fp,
				   const Projection* P,
				   const FormfactorPara* p,
				   double rmax, int npointsr,
				   int A, int l,
				   int jfin, int pfin, int afin,
				   int jini, int pini, int aini,
				   void* ffactorme,
				   const Eigenstates* E);

void writeTransitionChargeDensities(FILE* fp,
				    const Projection* P,
				    const FormfactorPara* p,
				    double rmax, int npointsr,
				    int Z, int l,
				    int jfin, int pfin, int afin,
				    int jini, int pini, int aini,
				    void* ffactorme,
				    const Eigenstates* E);

void writeTransitionChargeFourierBessel(FILE* fp,
					const Projection* P,
					const FormfactorPara* p,
					double rmax, int npointsr,
					int Z, int l,
					int jfin, int pfin, int afin,
					int jini, int pini, int aini,
					void* ffactorme,
					const Eigenstates* E);

#endif
