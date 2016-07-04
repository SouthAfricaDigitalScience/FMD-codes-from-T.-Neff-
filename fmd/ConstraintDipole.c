/**

  \file ConstraintDipole.c

  Constrain Dipole Moment


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
#include "ConstraintDipole.h"

#include "numerics/cmath.h"


Constraint ConstraintED2 = {
  val : 0.0,
  label : "sD2",
  me : calcConstraintED2,
  gradme : calcgradConstraintED2,
  output : outputConstraintED2
};


typedef struct {
  double Xcm[3];
  double d[3];
} dipolepara;


static void ob_edipole(const double* Xcm,
		       const Gaussian* G1, const Gaussian* G2, 
		       const GaussianAux* X, 
		       complex double d[3])
{
  int i;

  if (G1->xi == 1)
    for (i=0; i<3; i++)
      d[i] += (X->rho[i]-Xcm[i])* X->Q;
}


static void gob_edipole(const dipolepara* d2par,
			const Gaussian* G1, const Gaussian* G2,
			const GaussianAux* X, const gradGaussianAux* dX,
			complex double* ed2, gradGaussian* ded2)
{
  double* Xcm = d2par->Xcm;
  double* d = d2par->d;
  int i, k;
  complex double med;

  if (G1->xi == 1)
    for (i=0; i<3; i++) {
      med = X->rho[i]-Xcm[i];

      *ed2 += d[i]* med* X->Q;

      for (k=0; k<2; k++) 
	ded2->chi[k] += d[i]* med* dX->dQ.chi[k];
      ded2->a += d[i]*(dX->drho.a[i]* X->Q + med* dX->dQ.a);
      ded2->b[i] += d[i]* dX->drho.b* X->Q;
      for (k=0; k<3; k++)
	ded2->b[k] += d[i]* med* dX->dQ.b[k];
    }
}


static void gcmob_edipole(const dipolepara* d2par,
			const Gaussian* G1, const Gaussian* G2,
			const GaussianAux* X,
			gradCM* dCMed2)
{
  double* d = d2par->d;
  int i;

  if (G1->xi == 1)
    for (i=0; i<3; i++)
      dCMed2->X[i] -= d[i]* X->Q;
}


#define SQR(x) (x)*(x)

void calcConstraintED2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* ed2)
{
  double Xcm[3];
  double ed[3], edsq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_d = {dim: 3, opt: 1, par: Xcm, me: ob_edipole};

  calcSlaterDetOBME(Q, X, &op_ob_d, ed);

  for (i=0; i<3; i++)
    edsq += SQR(ed[i]);

  *ed2 = sqrt(edsq);
}


void calcgradConstraintED2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* ded2)
{
  double Xcm[3];
  double ed[3], edsq = 0.0;
  int i;
  dipolepara ed2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_d = {dim: 3, opt: 1, par: Xcm, me: ob_edipole};

  calcSlaterDetOBME(Q, X, &op_ob_d, ed);

  for (i=0; i<3; i++)
    edsq += SQR(ed[i]);
  
  for (i=0; i<3; i++)
    ed2par.Xcm[i] = Xcm[i];

  for (i=0; i<3; i++)
    ed2par.d[i] = ed[i]/sqrt(edsq);

  gradOneBodyOperator gop_ob_d2 = {opt: 1, par: &ed2par, me: gob_edipole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_d2, ded2);

  OneBodyOperator gcmop_ob_d2 = {dim: 6, opt: 1, par: &ed2par, me: gcmob_edipole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_d2, ded2);
  ded2->val = sqrt(edsq);
}


double outputConstraintED2(double val)
{
  return val;
}
