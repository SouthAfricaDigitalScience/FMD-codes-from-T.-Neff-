/**

   \file HOBasis.h

   define single-particle Harmonic oscillator basis


   (c) 2006 Thomas Neff

*/


#ifndef _HOBASIS_H
#define _HOBASIS_H

#include <fmd/Gaussian.h>


// implemented for basis up to Nmax = HONMAX 

#define HONMAX 8

// dim(harmonic oscillator)

#define HODIMMAX (HONMAX+1)*(HONMAX+2)*(HONMAX+3)/6


// harmonic oscillator single-particle states

typedef struct {
  int xi;	///< isospin [+1|-1]
  int n;	///< n=0,1,...
  int lj;	///< lj=0,1,2,...
  int m;	///< m=-2j,...,2j
} HOSPState;


void initHOBasis(int nmax);

int dimHOBasis(int nmax);

HOSPState* getHOspstate(int idx);
int HOspstateidx(const HOSPState* spstate);
int HOxinljmidx(int xi, int n, int l, int twoj, int twom);

char* labelHOshell(int N);
char* labelHOorbit(int n, int l, int twoj);
char* labelHOspstate(const HOSPState* spstate);

// lj index into 2*l, 2*j

inline static int ljtotwol(int lj) { return 2*((lj+1)/2); }
inline static int ljtotwoj(int lj) { return 2*(lj/2)+1; }

// calculate the overlap of Gaussian with all HO single-particle states

void amplitudesHOGaussian(const Gaussian* G, 
			  int nmax, double omega,
			  complex double* amp); 


#endif
