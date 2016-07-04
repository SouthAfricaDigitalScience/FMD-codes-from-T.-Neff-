/**

  \file gradCMSlaterDet.h

  gradients of matrix elements that that depend
  on the expectation values of Center of Mass coordinate and momentum


  (c) 2003 Thomas Neff

*/


#ifndef _GRADCMSLATERDET_H
#define _GRADCMSLATERDET_H


#include "Gaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"


typedef struct {
  complex double X[3];
  complex double V[3];
} gradCM;


void calcgradCMSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const OneBodyOperator* op,
			     gradSlaterDet* grad);


void calcgradCMSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const TwoBodyOperator* op,
			     gradSlaterDet* grad);

#endif
