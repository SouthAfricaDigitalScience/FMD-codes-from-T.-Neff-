/**

  \file ConstraintOctupole.c

  Constrain Mass Octupole
  center of mass corrected by substraction of Xcm expectation value


  (c) 2003 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "CenterofMass.h"
#include "gradCMSlaterDet.h"

#include "Constraint.h"
#include "ConstraintOctupole.h"

#include "numerics/cmath.h"


Constraint ConstraintO2 = {
  val : 0.0,
  label : "sO2",
  me : calcConstraintO2,
  gradme : calcgradConstraintO2,
  output : outputConstraintO2
};

Constraint ConstraintEO2 = {
  val : 0.0,
  label : "sEO2",
  me : calcConstraintEO2,
  gradme : calcgradConstraintEO2,
  output : outputConstraintEO2
};

Constraint ConstraintNO2 = {
  val : 0.0,
  label : "sNO2",
  me : calcConstraintNO2,
  gradme : calcgradConstraintNO2,
  output : outputConstraintNO2
};


typedef struct {
  double Xcm[3];
  double o[27];
} octupolepara;


static void ob_octupole(const double* Xcm,
			const Gaussian*G1, const Gaussian* G2, 
			const GaussianAux* X, 
			complex double o[27])
{
  complex double rho2cm = 0.0;
  int i,j,k;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i] - Xcm[i]);

  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++) {
	o[i+j*3+k*9] += 
	  5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])*X->Q;

	if (i==j) o[i+j*3+k*9] += - rho2cm* (X->rho[k]-Xcm[k])* X->Q; 
	if (j==k) o[i+j*3+k*9] += - rho2cm* (X->rho[i]-Xcm[i])* X->Q; 
	if (k==i) o[i+j*3+k*9] += - rho2cm* (X->rho[j]-Xcm[j])* X->Q; 
      }
}


static void ob_eoctupole(const double* Xcm,
			 const Gaussian*G1, const Gaussian* G2, 
			 const GaussianAux* X, 
			 complex double o[27])
{
  complex double rho2cm = 0.0;
  int i,j,k;

  if (G1->xi == 1) {
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i] - Xcm[i]);

    for (k=0; k<3; k++)
      for (j=0; j<3; j++)
	for (i=0; i<3; i++) {
	  o[i+j*3+k*9] += 
	    5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])*X->Q;

	  if (i==j) o[i+j*3+k*9] += - rho2cm* (X->rho[k]-Xcm[k])* X->Q; 
	  if (j==k) o[i+j*3+k*9] += - rho2cm* (X->rho[i]-Xcm[i])* X->Q; 
	  if (k==i) o[i+j*3+k*9] += - rho2cm* (X->rho[j]-Xcm[j])* X->Q; 
	}
  }
}


static void ob_noctupole(const double* Xcm,
			 const Gaussian*G1, const Gaussian* G2, 
			 const GaussianAux* X, 
			 complex double o[27])
{
  complex double rho2cm = 0.0;
  int i,j,k;

  if (G1->xi == -1) {
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i] - Xcm[i]);

    for (k=0; k<3; k++)
      for (j=0; j<3; j++)
	for (i=0; i<3; i++) {
	  o[i+j*3+k*9] += 
	    5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])*X->Q;

	  if (i==j) o[i+j*3+k*9] += - rho2cm* (X->rho[k]-Xcm[k])* X->Q; 
	  if (j==k) o[i+j*3+k*9] += - rho2cm* (X->rho[i]-Xcm[i])* X->Q; 
	  if (k==i) o[i+j*3+k*9] += - rho2cm* (X->rho[j]-Xcm[j])* X->Q; 
	}
  }
}


static void gob_octupole(const octupolepara* P,
			 const Gaussian* G1, const Gaussian* G2,
			 const GaussianAux* X, const gradGaussianAux* dX,
			 complex double* o2, gradGaussian* do2)
{
  double* Xcm = P->Xcm;
  double* o=P->o;
  int i,j,k,l;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;
  complex double merrr;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i] - Xcm[i]);

  drho2cm.a = 0.0;
  for (i=0; i<3; i++)
    drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
  for (i=0; i<3; i++)
    drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++) {

	merrr = 5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]);
      
	*o2 += o[i+j*3+k*9]* merrr* X->Q;

	for (l=0; l<2; l++)
	  do2->chi[l] +=  o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	do2->a += o[i+j*3+k*9]*
	  (5*(dX->drho.a[i]*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]) + 
	      (X->rho[i]-Xcm[i])*dX->drho.a[j]*(X->rho[k]-Xcm[k]) +
	      (X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.a[k])* X->Q+
	   merrr* dX->dQ.a);
	do2->b[i] += o[i+j*3+k*9]*
	  5*dX->drho.b*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])* X->Q;
	do2->b[j] += o[i+j*3+k*9]*
	  5*(X->rho[i]-Xcm[i])*dX->drho.b*(X->rho[k]-Xcm[k])* X->Q;
	do2->b[k] += o[i+j*3+k*9]*
	  5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.b* X->Q;
	for (l=0; l<3; l++) 
	  do2->b[l] += o[i+j*3+k*9]* merrr* dX->dQ.b[l];
 
	if (i==j) {
	  merrr = rho2cm*(X->rho[k]-Xcm[k]);
	  *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	  for (l=0; l<2; l++)
	    do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	  do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[k]-Xcm[k])*X->Q + 
				    rho2cm*dX->drho.a[k]*X->Q +
				    merrr* dX->dQ.a);
	  do2->b[k] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	  for (l=0; l<3; l++)
	    do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[k]-Xcm[k])*X->Q+
					 merrr* dX->dQ.b[l]);
	}

	if (i==k) {
	  merrr = rho2cm*(X->rho[j]-Xcm[j]);
	  *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	  for (l=0; l<2; l++)
	    do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	  do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[j]-Xcm[j])*X->Q + 
				    rho2cm*dX->drho.a[j]*X->Q +
				    merrr* dX->dQ.a);
	  do2->b[j] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	  for (l=0; l<3; l++)
	    do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[j]-Xcm[j])*X->Q+
					 merrr* dX->dQ.b[l]);
	}

	if (j==k) {
	  merrr = rho2cm*(X->rho[i]-Xcm[i]);
	  *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	  for (l=0; l<2; l++)
	    do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	  do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[i]-Xcm[i])*X->Q + 
				    rho2cm*dX->drho.a[i]*X->Q +
				    merrr* dX->dQ.a);
	  do2->b[i] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	  for (l=0; l<3; l++)
	    do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[i]-Xcm[i])*X->Q+
					 merrr* dX->dQ.b[l]);
	}

  }
}


static void gob_eoctupole(const octupolepara* P,
			  const Gaussian* G1, const Gaussian* G2,
			  const GaussianAux* X, const gradGaussianAux* dX,
			  complex double* o2, gradGaussian* do2)
{
  double* Xcm = P->Xcm;
  double* o=P->o;
  int i,j,k,l;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;
  complex double merrr;

  if (G1->xi == 1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i] - Xcm[i]);

    drho2cm.a = 0.0;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    for (k=0; k<3; k++)
      for (j=0; j<3; j++)
	for (i=0; i<3; i++) {

	  merrr = 5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]);
      
	  *o2 += o[i+j*3+k*9]* merrr* X->Q;

	  for (l=0; l<2; l++)
	    do2->chi[l] +=  o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	  do2->a += o[i+j*3+k*9]*
	    (5*(dX->drho.a[i]*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]) + 
		(X->rho[i]-Xcm[i])*dX->drho.a[j]*(X->rho[k]-Xcm[k]) +
		(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.a[k])* X->Q+
	     merrr* dX->dQ.a);
	  do2->b[i] += o[i+j*3+k*9]*
	    5*dX->drho.b*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])* X->Q;
	  do2->b[j] += o[i+j*3+k*9]*
	    5*(X->rho[i]-Xcm[i])*dX->drho.b*(X->rho[k]-Xcm[k])* X->Q;
	  do2->b[k] += o[i+j*3+k*9]*
	    5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.b* X->Q;
	  for (l=0; l<3; l++) 
	    do2->b[l] += o[i+j*3+k*9]* merrr* dX->dQ.b[l];
 
	  if (i==j) {
	    merrr = rho2cm*(X->rho[k]-Xcm[k]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[k]-Xcm[k])*X->Q + 
				      rho2cm*dX->drho.a[k]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[k] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[k]-Xcm[k])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	  if (i==k) {
	    merrr = rho2cm*(X->rho[j]-Xcm[j]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[j]-Xcm[j])*X->Q + 
				      rho2cm*dX->drho.a[j]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[j] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[j]-Xcm[j])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	  if (j==k) {
	    merrr = rho2cm*(X->rho[i]-Xcm[i]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[i]-Xcm[i])*X->Q + 
				      rho2cm*dX->drho.a[i]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[i] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[i]-Xcm[i])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	}
  }
}


static void gob_noctupole(const octupolepara* P,
			  const Gaussian* G1, const Gaussian* G2,
			  const GaussianAux* X, const gradGaussianAux* dX,
			  complex double* o2, gradGaussian* do2)
{
  double* Xcm = P->Xcm;
  double* o=P->o;
  int i,j,k,l;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;
  complex double merrr;

  if (G1->xi == -1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i] - Xcm[i]);

    drho2cm.a = 0.0;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    for (k=0; k<3; k++)
      for (j=0; j<3; j++)
	for (i=0; i<3; i++) {

	  merrr = 5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]);
      
	  *o2 += o[i+j*3+k*9]* merrr* X->Q;

	  for (l=0; l<2; l++)
	    do2->chi[l] +=  o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	  do2->a += o[i+j*3+k*9]*
	    (5*(dX->drho.a[i]*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k]) + 
		(X->rho[i]-Xcm[i])*dX->drho.a[j]*(X->rho[k]-Xcm[k]) +
		(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.a[k])* X->Q+
	     merrr* dX->dQ.a);
	  do2->b[i] += o[i+j*3+k*9]*
	    5*dX->drho.b*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])* X->Q;
	  do2->b[j] += o[i+j*3+k*9]*
	    5*(X->rho[i]-Xcm[i])*dX->drho.b*(X->rho[k]-Xcm[k])* X->Q;
	  do2->b[k] += o[i+j*3+k*9]*
	    5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*dX->drho.b* X->Q;
	  for (l=0; l<3; l++) 
	    do2->b[l] += o[i+j*3+k*9]* merrr* dX->dQ.b[l];
 
	  if (i==j) {
	    merrr = rho2cm*(X->rho[k]-Xcm[k]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[k]-Xcm[k])*X->Q + 
				      rho2cm*dX->drho.a[k]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[k] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[k]-Xcm[k])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	  if (i==k) {
	    merrr = rho2cm*(X->rho[j]-Xcm[j]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[j]-Xcm[j])*X->Q + 
				      rho2cm*dX->drho.a[j]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[j] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[j]-Xcm[j])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	  if (j==k) {
	    merrr = rho2cm*(X->rho[i]-Xcm[i]);
	    *o2 += - o[i+j*3+k*9]* merrr* X->Q;

	    for (l=0; l<2; l++)
	      do2->chi[l] += - o[i+j*3+k*9]* merrr* dX->dQ.chi[l];
	    do2->a += - o[i+j*3+k*9]*(drho2cm.a*(X->rho[i]-Xcm[i])*X->Q + 
				      rho2cm*dX->drho.a[i]*X->Q +
				      merrr* dX->dQ.a);
	    do2->b[i] += -o[i+j*3+k*9]* rho2cm* dX->drho.b*X->Q;
	    for (l=0; l<3; l++)
	      do2->b[l] += - o[i+j*3+k*9]*(drho2cm.b[l]*(X->rho[i]-Xcm[i])*X->Q+
					   merrr* dX->dQ.b[l]);
	  }

	}
  }
}


static void gcmob_octupole(const octupolepara* P,
			   const Gaussian* G1, const Gaussian* G2, 
			   const GaussianAux* X, 
			   gradCM* dCMo2)
{
  double* Xcm = P->Xcm;
  double* o = P->o;
  complex double rho2cm = 0.0;
  gradCM drho2cm;
  int i,j,k,l;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  for (i=0; i<3; i++)
    drho2cm.X[i] = -2*(X->rho[i]-Xcm[i]);

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      for (k=0; k<3; k++) {
	dCMo2->X[i] -= o[i+j*3+k*9]*
	  5*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])*X->Q;
	dCMo2->X[j] -= o[i+j*3+k*9]*
	  5*(X->rho[i]-Xcm[i])*(X->rho[k]-Xcm[k])*X->Q;
	dCMo2->X[k] -= o[i+j*3+k*9]*
	  5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;

	if (i==j) {
	  for (l=0; l<3; l++)
	    dCMo2->X[l] += - o[i+j*3+k*9]*drho2cm.X[l]*(X->rho[k]-Xcm[k])* X->Q;
	  dCMo2->X[k] += o[i+j*3+k*9]*rho2cm* X->Q;
	}
	if (j==k) {
	  for (l=0; l<3; l++)
	    dCMo2->X[l] += - o[i+j*3+k*9]*drho2cm.X[l]*(X->rho[i]-Xcm[i])* X->Q;
	  dCMo2->X[i] += o[i+j*3+k*9]*rho2cm* X->Q;
	}
	if (k==i) {
	  for (l=0; l<3; l++)
	    dCMo2->X[l] += - o[i+j*3+k*9]*drho2cm.X[l]*(X->rho[j]-Xcm[j])* X->Q;
	  dCMo2->X[j] += o[i+j*3+k*9]*rho2cm* X->Q;
	}
      }
}


#define SQR(x) (x)*(x)

void calcConstraintO2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* o2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_octupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);

  *o2 = sqrt(osq);
}


void calcConstraintEO2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* o2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_eoctupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);

  *o2 = sqrt(osq);
}


void calcConstraintNO2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* o2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_noctupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);

  *o2 = sqrt(osq);
}


void calcgradConstraintO2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* do2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;
  octupolepara o2para;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_octupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);
  
  for (i=0; i<3; i++)
    o2para.Xcm[i] = Xcm[i];

  for (i=0; i<27; i++)
    o2para.o[i] = o[i]/sqrt(osq);

  gradOneBodyOperator gop_ob_o2 = {opt: 1, par: &o2para, me: gob_octupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_o2, do2);

  OneBodyOperator gcmop_ob_o2 = {dim: 6, opt: 1, par: &o2para, me: gcmob_octupole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_o2, do2);
  do2->val = sqrt(osq);
}

void calcgradConstraintEO2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* do2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;
  octupolepara o2para;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_eoctupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);
  
  for (i=0; i<3; i++)
    o2para.Xcm[i] = Xcm[i];

  for (i=0; i<27; i++)
    o2para.o[i] = o[i]/sqrt(osq);

  gradOneBodyOperator gop_ob_o2 = {opt: 1, par: &o2para, me: gob_eoctupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_o2, do2);

  OneBodyOperator gcmop_ob_o2 = {dim: 6, opt: 1, par: &o2para, me: gcmob_octupole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_o2, do2);
  do2->val = sqrt(osq);
}

void calcgradConstraintNO2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* do2)
{
  double Xcm[3];
  double o[27], osq = 0.0;
  int i;
  octupolepara o2para;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_o = {dim: 27, opt: 1, par: Xcm, me: ob_noctupole};

  calcSlaterDetOBME(Q, X, &op_ob_o, o);

  for (i=0; i<27; i++)
    osq += SQR(o[i]);
  
  for (i=0; i<3; i++)
    o2para.Xcm[i] = Xcm[i];

  for (i=0; i<27; i++)
    o2para.o[i] = o[i]/sqrt(osq);

  gradOneBodyOperator gop_ob_o2 = {opt: 1, par: &o2para, me: gob_noctupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_o2, do2);

  OneBodyOperator gcmop_ob_o2 = {dim: 6, opt: 1, par: &o2para, me: gcmob_octupole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_o2, do2);
  do2->val = sqrt(osq);
}


double outputConstraintO2(double val)
{
  return val;
}


double outputConstraintEO2(double val)
{
  return val;
}


double outputConstraintNO2(double val)
{
  return val;
}


