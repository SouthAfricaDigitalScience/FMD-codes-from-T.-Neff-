/*

  cmat.h

  complex matrices


  (c) 2003,2004 Thomas Neff

*/


#ifndef _CMAT_H
#define _CMAT_H


#include <stdio.h>
#include <complex.h>


/// initialize complex nxn matrix
complex double* initcmat(int n);


/// copy matrix A into B
void copycmat(int n, const complex double* A, complex double* B);


/// mirror upper triangle of hermitian matrix into lower triangle
void mirrorcmat(int n, complex double* A);


/// make matrix hermitian
void hermitizecmat(int n, complex double* A);


/// multiply complex matrices C = A.B
void multcmat(const complex double* A, const complex double* B,
	      complex double* C, int n);


/// calculate pseudoinverse B of A using SVD
void pseudoinverse(const complex double* A, complex double* B,
		   int n, double thresh);


/// solve eigenvalue problem for complex matrix A
void eigensystem(complex double* A,
		 complex double* a, complex double* V, int n);

/// solve generalized eigenvalue problem for complex matrices A and N
/// using SVD
void generalizedeigensystem(complex double* A, complex double* N, int n,
			    double thresh,
			    complex double* v, complex double* V, int* dim);

/// sort first dim eigenstates, starting with smallest eigenvalue
void sorteigenstates(int n, complex double* v, complex double* V, int dim);

/// read integer vector from fp
void freadivec(FILE* fp, int n, int* v);

/// write vector to fp
void fprintcvec(FILE* fp, int n, const complex double* V);

/// write vector to fp in binary format
void fwritecvecbin(FILE* fp, int n, const complex double* V);

/// read vector from fp
void freadcvec(FILE* fp, int n, complex double* V);


/// write matrix to fp
void fprintcmat(FILE* fp, int n, const complex double* A);

/// write matrix to fp in binary format
void fwritecmatbin(FILE* fp, int n, const complex double* A);

/// read matrix from fp
void freadcmat(FILE* fp, int n, complex double* A);


/// write colums of matrix to fp
void fprintcmatcols(FILE* fp, int n, int cols, const complex double* A);

/// read colums of matrix from fp
void freadcmatcols(FILE* fp, int n, int cols, complex double* A);


#endif
