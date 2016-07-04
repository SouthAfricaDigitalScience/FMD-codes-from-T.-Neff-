/**

  \file Projection.h

  angular momentum and linear momentum projection of matrix elements


  (c) 2005-2011 Thomas Neff

*/


#ifndef _PROJECTION_H
#define _PROJECTION_H

#include <zlib.h>

#include "SlaterDet.h"
#include "Interaction.h"
#include "Observables.h"
#include "Symmetry.h"


/// Integration for CM projection using Product integration or Polyhedron points
enum { CMNone, CMSimple, CMProd, CMPoly };

typedef struct {
  unsigned int nr     : 8;
  unsigned int ntheta : 8;
  unsigned int nphi   : 8;
} cmprodpara;

typedef struct {
  unsigned int nr     : 8;
  unsigned int npoly  : 8;
} cmpolypara;


typedef struct {
  int type;
  int n;
  union {
    cmprodpara prod;
    cmpolypara poly;
  };
  double (*x)[3];
  double* w;
} cmintegrationpara;


/// Integration over Euler angles using Product integration or ZCW points
enum { AngNone, AngProd, AngProdA, AngZCW };

typedef struct {
  unsigned int nazimuth : 8;
  unsigned int nbeta  : 8;
} angprodpara;

typedef struct {
  unsigned int idx : 8;
} angzcwpara;


typedef struct {
  int idx;
} angzcwintegrationpara;

typedef struct {
  int nalpha;
  double* alpha;
  double* walpha;
  int nbeta;
  double* beta;
  double* wbeta;
  int ngamma;
  double* gamma;
  double* wgamma;
} angprodintegrationpara;


typedef struct {
  int type;
  int n;
  union {
    angzcwintegrationpara zcw;
    angprodintegrationpara prod;
  };   
} angintegrationpara;


/// Projection
typedef struct {
  int jmax;		///< calculate up to (jmax-1)/2
  int odd;		///< odd particle number - half integer spin
  //  this is the property of right state in the
  //  matrix element, oddness of left state has to be deduced !
  int ang;		///< type of angular momentum projection
  union {
    angprodpara angprod;
    angzcwpara angzcw;
  };
  int cm;		///< type of cm projection
  union {
    cmprodpara cmprod;
    cmpolypara cmpoly;
  };
} Projection;


/// General ManyBody Operator
typedef struct {
  char* name;		///< uniquely identify operator including parameters
  int rank;		///< tensor rank of operator, 0 scalar, 2 vector, ...
  int pi;		///< parity of operator, 0 positve, 1 negative
  int dim;		///< dimension of operator
  int size;		///< might be bigger than dimension
  void* par;
  void (*me)(void* parameters,
	     const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X,
	     complex double *val);
} ManyBodyOperator;


/// Collection of ManyBody Operators
typedef struct {
  char* name;
  int n;		///< number of operators
  void* par;		///< the same for all operators
  int dim;		///< has to be the same for all Operators
  int size;
  ManyBodyOperator* Op;
  void (*me)(void* parameters,
	     const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X,
	     complex double **val);
} ManyBodyOperators;


/// Container for Eigenstates
typedef struct {
  int n;
  complex double** v;		///< energy eigenvalues
  complex double** V;		///< eigenvectors
  int* dim;			///< number of eigenstates calculated by SVD
  int* ngood;			///< number of eigenstates considered as "good"
  int** index;			///< indexing eigenstates according to "goodness"
  complex double** norm;	///< many-body norm^2 of eigenstates
} Eigenstates;


typedef struct {
  int n;
  int** ngood;
  complex double** amp;
} Amplitudes;


// default jmax value
#define JMAX 9

// integer j is running 0, 2, 4, 6, ...
// half integer j is running 1, 3, 5, 7, ...

inline static int idxpij(int jmax, int pi, int j)
{
  return(j/2 + pi*(jmax+1)/2);
}

// integer j is running 0, 2, 4, 6, 8
// half integer j is running 1, 3, 5, 7

inline static int idxjm(int j, int m)
{
  return ((j+m)/2);
}


inline static int idxnjm(int a, int j, int m)
{
  return (idxjm(j,m)+a*(j+1));
}


inline static int idxjmk(int j, int m, int k)
{
  return (idxjm(j,m) + idxjm(j,k)*(j+1));
}

char* AngmomtoStr(int j, int pi);

char* ProjectiontoStr(const Projection* P);

int initAngintegration(angintegrationpara* par, const char* projpar);

void freeAngintegration(angintegrationpara* par);

int initProjection(Projection* P, int odd, const char* projpar);


void fprintProjectinfo(FILE* fp, const Projection* P);


void* initprojectedVector(const Projection* P, const ManyBodyOperator* Op, 
			  int n);

void* initprojectedVectornull(const Projection* P, const ManyBodyOperator* Op, 
			      int n);


void* initprojectedtransitionVector(const Projection* P, const ManyBodyOperator* Op,
				    int n, int np);

void* initprojectedtransitionVectornull(const Projection* P, const ManyBodyOperator* Op,
					int n, int np);


void* initprojectedMBME(const Projection* P, const ManyBodyOperator* Op);


void initcmintegration(const Projection* P, 
		       const SlaterDet* Q, const SlaterDet* Qp,
		       cmintegrationpara* par);

double _estimateacm(const SlaterDet* Q);

void initCMintegration(cmintegrationpara* par, const char* projpar);

void _initcmintegration(const Projection* P,
                        double alpha,
                        cmintegrationpara* par);

void reinitcmintegration(const SlaterDet* Q, const SlaterDet* Qp,
			 cmintegrationpara* par);

void freecmintegration(cmintegrationpara* par);

void getcmintegrationpoint(int i, const cmintegrationpara* par,
			   double x[3], double* w);


void initangintegration(const Projection* P, 
			const SlaterDet* Q, const SlaterDet* Qp,
			Symmetry S, Symmetry Sp,
			angintegrationpara* par);

double _estimateangkappa(const SlaterDet* Q);

void _setangkappacrit(double kappa);

double _getangkappacrit();

void _initangintegration(const Projection* P,
                         double kappa,
                         Symmetry S, Symmetry Sp,
                         angintegrationpara* par);


void reinitangintegration(const SlaterDet* Q, const SlaterDet* Qp,
                          Symmetry S, Symmetry Sp,
                          angintegrationpara* par);


void getangintegrationpoint(int i, const angintegrationpara* par,
			    double* alpha, double* beta, double* gamma, 
			    double* w);


void calcprojectedMBME(const Projection* P, const ManyBodyOperator* Op,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       Symmetry S, Symmetry Sp,
		       void* mbme);


void calcprojectedMBMEs(const Projection* P, const ManyBodyOperators* Ops,
			const SlaterDet* Q, const SlaterDet* Qp,
			Symmetry S, Symmetry Sp,
			void* mbmes);


void calcprojectedMBMEsingular(const Projection* P, const ManyBodyOperator* Op,
			       const SlaterDet* Q, const SlaterDet* Qp,
			       Symmetry S, Symmetry Sp,
			       void* mbme);


void hermitizeprojectedMBME(const Projection* P, const ManyBodyOperator* Op,
			    void* mbme, int n);


int writeprojectedMBME(gzFile fp,
		       const Projection* P, const ManyBodyOperator* Op,
		       Symmetry S, Symmetry Sp,
		       const void* mbme);


int writeprojectedMBMEtoFile(const char* mbfile, const char* mbfilep,
			     const Projection* P, const ManyBodyOperator* Op,
			     Symmetry S, Symmetry Sp,
			     const void* mbme);


int readprojectedMBME(gzFile fp,
		      const Projection* P, const ManyBodyOperator* Op,
		      Symmetry S, Symmetry Sp,
		      void* mbme);


int readprojectedMBMEfromFile(const char* mbfile, const char* mbfilep,
			      const Projection* P, const ManyBodyOperator* Op,
			      Symmetry S, Symmetry Sp,
			      void* mbme);


void initEigenstates(const Projection* P, Eigenstates* E, int n);


void initAmplitudes(const Projection* P, Amplitudes* A, int n);


void calcEigenstates(const Projection* P, 
		     const Interaction* Int,
		     const Observablesod*** obsme,
		     Eigenstates* E, 
		     double thresh);


void calcEigenstatesK(const Projection* P, 
		      const Interaction* Int,
		      const Observablesod*** obsme,
		      Eigenstates* E, 
		      int K, double thresh);



void calcMultiEigenstates(const Projection* P,
			  const Interaction* Int,
			  const Observablesod*** obsme,
			  const Eigenstates* Ep,
			  Eigenstates* multiE, Amplitudes* multiA, 
			  double thresh);


/// calculate expectation values with matrix elements mbme
/// and eigenvectors E
void calcexpectprojectedMBME(const Projection* P,
			     const ManyBodyOperator* Op,
			     const void* mbme,
			     const Symmetry* S,
			     const Eigenstates* E,
			     void* expectmbme);

void calcexpectprojectedMBMEipj(const Projection* P,
				const ManyBodyOperator* Op,
				const void* mbme,
				const Symmetry* S,
				const Eigenstates* E,
				int j, int p, int i,
				void* expectmbme);


/// calculate transition strengths with matrix elements mbme
/// and eigenvectors E
void calctransitionprojectedMBME(const Projection* P,
				 const ManyBodyOperator* Op,
				 const void* mbme,
				 const Symmetry* S,
				 const Symmetry* Sp,
				 const Eigenstates* E,
				 const Eigenstates* Ep,
				 void* transmbme);

void calctransitionprojectedMBMEipj(const Projection* P,
				    const ManyBodyOperator* Op,
				    const void* mbme,
				    const Symmetry* S,
				    const Symmetry* Sp,
				    const Eigenstates* E,
				    const Eigenstates* Ep,
				    int jini, int pini, int iini,
				    int jfin, int pfin, int ifin,
				    void* transmbme);


void sortEigenstates(const Projection* P,
		     const Interaction* Int,
		     const Observablesod** obs,
		     Eigenstates* E, 
		     double minnorm, int all);



int writeEigenstates(FILE* fp, const Projection* P, const Eigenstates* E);


int readEigenstates(FILE* fp, const Projection* P, Eigenstates* E, int n);


int readEigenstatesfromFile(const char* fname,
			    const Projection* P, Eigenstates* E, int n);


int writeMulticonfigfile(const char* fname,
			 const Projection* P,
			 Symmetry* S,
			 const char** mbfile,
			 const Eigenstates* E,
			 int n);


int readMulticonfigfile(const char* fname,
			char*** slaterdetfile,
			Projection* P,
			SlaterDet** Q,
			Symmetry** S,
			Eigenstates* E,
			int* nstates);

#endif
