/**

  \file ConstraintSymmetries.h

  constrain intrinsic state regarding discrete symmetries


  (c) 2009 Thomas Neff

*/


#ifndef _CONSTRAINTSYMMETRIES_H
#define _CONSTRAINTSYMMETRIES_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Constraint.h"


extern Constraint ConstraintParityP;

extern Constraint ConstraintParityN;

extern Constraint ConstraintParityCharge;

extern Constraint ConstraintD3H;

extern Constraint ConstraintSxz;

extern Constraint ConstraintSyz;


void calcConstraintParP(const SlaterDet* Q, const SlaterDetAux* X,
                        double* parpconstraint);

void calcgradConstraintParP(const SlaterDet* Q, const SlaterDetAux* X,
                            const gradSlaterDetAux* dX,
                            gradSlaterDet* dparpconstraint);


void calcConstraintParN(const SlaterDet* Q, const SlaterDetAux* X,
                        double* parnconstraint);

void calcgradConstraintParN(const SlaterDet* Q, const SlaterDetAux* X,
                            const gradSlaterDetAux* dX,
                            gradSlaterDet* dparnconstraint);


void calcConstraintParCh(const SlaterDet* Q, const SlaterDetAux* X,
                         double* parchconstraint);

void calcgradConstraintParCh(const SlaterDet* Q, const SlaterDetAux* X,
                             const gradSlaterDetAux* dX,
                             gradSlaterDet* dparchconstraint);


void calcConstraintD3H(const SlaterDet* Q, const SlaterDetAux* X,
                       double* d3hconstraint);

void calcgradConstraintD3H(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dd3hconstraint);


void calcConstraintSxz(const SlaterDet* Q, const SlaterDetAux* X,
                       double* sxzconstraint);

void calcgradConstraintSxz(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dsxzconstraint);


void calcConstraintSyz(const SlaterDet* Q, const SlaterDetAux* X,
                       double* syzconstraint);

void calcgradConstraintSyz(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dsyzconstraint);


double outputConstraintZero(double val);


#endif
