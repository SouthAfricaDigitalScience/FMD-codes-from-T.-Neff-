/**

  \file Constraint.h

  defining a constraint


  (c) 2003 Thomas Neff

*/


#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include "SlaterDet.h"
#include "gradSlaterDet.h"


typedef struct {
  double val;
  char* label;
  void (*me)(const SlaterDet* Q, const SlaterDetAux* X, double* val);
  void (*gradme)(const SlaterDet* Q, const SlaterDetAux* X,
		 const gradSlaterDetAux* dX, gradSlaterDet* dval);
  double (*output)(double val);
} Constraint;


typedef struct {
  double val;
  char* label;
  int dim;
  int gradneedsme;
  void (*me)(const SlaterDet* Q, const SlaterDet* Qp, 
	     const SlaterDetAux* X, complex double* mes);
  void (*gradme)(const SlaterDet* Q, const SlaterDet* Qp, 
		 const SlaterDetAux* X,
		 const gradSlaterDetAux* dX, 
		 const double* mes, 
		 gradSlaterDet* dval);
  double (*output)(double val);
} Constraintod;


#endif
