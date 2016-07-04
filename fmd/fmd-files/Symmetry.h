/**

  \file Symmetry.h

  Symmetries and Selection Rules for intrinsic states


  (c) 2005 Thomas Neff

*/


#ifndef _SYMMETRY_H
#define _SYMMETRY_H


/// encode the symmetries as bits in an integer

enum { parity0=0, parity1, spherical,
       axial=3,					  // indicates any of axialx
       axial0=4, axial1, axial2, axial3, axial4, 
       axial5, axial6, axial7, axial8, axial9,
       rotatez2=16, rotatez3, 
       rotatey1=20, 
       reflectxy=24,
       reflectxz };

typedef unsigned int Symmetry;


inline static int hasSymmetry(Symmetry S, int s)
{
  return (S >> s & 1);
}

inline static void setSymmetry(Symmetry* S, int s)
{
  *S |= 1 << s;
}


// not very elegant, should use bitwise and
inline static int hasAxialSymmetry(Symmetry S)
{
  return (hasSymmetry(S, spherical) ||
	  hasSymmetry(S, axial)  ||
	  hasSymmetry(S, axial0) ||
	  hasSymmetry(S, axial1) ||
	  hasSymmetry(S, axial2) ||
	  hasSymmetry(S, axial3) ||
	  hasSymmetry(S, axial4) ||
	  hasSymmetry(S, axial5) ||
	  hasSymmetry(S, axial6) ||
	  hasSymmetry(S, axial7) ||
	  hasSymmetry(S, axial8) ||
	  hasSymmetry(S, axial9));
}

/// quantum numbers J^pi,K allowed by symmetry ?
int SymmetryAllowed(Symmetry S, int pi, int j, int k);

/// get Symmetry from string 'SYMMETRY:FILENAME'
/// return pointer to FILENAME
void extractSymmetryfromString(char** str, Symmetry* S);

/// get String with Symmetry
char* SymmetrytoStr(Symmetry S);


#endif
