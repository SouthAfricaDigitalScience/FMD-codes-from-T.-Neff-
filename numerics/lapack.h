/*

  lapack.h

  declarations for Lapack routines


  (c) 2003 Thomas Neff

*/


#ifndef _LAPACK_H
#define _LAPACK_H

#include <complex.h>
#include "fortranc.h"


/// compute the solution to a real system of linear equations A*X = B
void FORTRAN(dposv)(const char* UPLO, const int* N, const int* NRHS,
		    double* A, const int* LDA, double* B, const int* LDB,
		    int* INFO);


/// compute all eigenvalues and, optionally, eigenvectors of a real
/// symmetric matrix A
void FORTRAN(dsyev)(const char* JOBZ, const char* UPLO, const int* N,
		    double* A, const int* LDA, double* W, double* WORK,
		    const int* LWORK, int* INFO);


/// cholesky factorization of complex Hermitian positive definite matrix A
void FORTRAN(zpotrf)(const char* UPLO, const int* N, complex double* A,
		     const int* LDA, int* INFO);


/// inverse of complex Hermitian positive definite matrix A using Cholesky
/// factorization computed by zpotrf
void FORTRAN(zpotri)(const char* UPLO, const int* N, complex double* A,
		     const int* LDA, int* INFO);


/// compute an LU factorization of a general M-by-N matrix A
void FORTRAN(zgetrf)(const int* M, const int* N, complex double* A, 
		     const int* LDA, int* IPIV, int* INFO);


/// computes the inverse of a matrix using the LU factorization 
/// computed  by  zgetrf 
void FORTRAN(zgetri)(const int* N, complex double* A, 
		     const int* LDA, int* IPIV, 
		     complex double* WORK, const int* LWORK, int* INFO);


/// calculate determinant of matrix using the LU factorization 
/// computed  by  zgetrf 
void FORTRAN(zdet)(const complex double* A, const int* LDA, 
		   const int* N, int* IPIV, complex double* DET);


/// compute the singular value decomposition (SVD) of a complex matrix
void FORTRAN(zgesvd)(const char* JOBU, const char* JOBVT, 
		     const int* M, const int* N,
		     complex double* A, const int* LDA, double* S,
		     complex double* U, const int* LDU, 
		     complex double* VT, const int* LDVT, 
		     complex double* WORK, const int* LWORK, 
		     double* RWORK, int* INFO);


/// compute the eigenvalues and eigenvectors of complex matrix A
void FORTRAN(zgeev)(const char* JOBVL, const char* JOBVR, 
		    const int* N, complex double* A, const int* LDA, 
		    complex double* W, complex double* VL, 
		    const int* LDVL, complex double* VR, const int* LDVR,
		    complex double* WORK, const int* LWORK, double* RWORK, 
		    int* INFO);

#endif
