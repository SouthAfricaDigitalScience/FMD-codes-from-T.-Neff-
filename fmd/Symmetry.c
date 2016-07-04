/**

  \file Symmetry.c

  Symmetries and Selection Rules for intrinsic states


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <string.h>


#include "Symmetry.h"


int SymmetryAllowed(Symmetry sym, int pi, int j, int k)
{
  if (sym==0) return 1;

  if (hasSymmetry(sym, spherical) && (j!=0 || pi==1)) return 0;
  if (hasSymmetry(sym, parity0) && pi==1) return 0;
  if (hasSymmetry(sym, parity1) && pi==0) return 0;
  if (hasSymmetry(sym, axial0) && k!=0) return 0;
  if (hasSymmetry(sym, axial1) && k!=1) return 0;
  if (hasSymmetry(sym, axial2) && k!=2) return 0;
  if (hasSymmetry(sym, axial3) && k!=3) return 0;
  if (hasSymmetry(sym, axial4) && k!=4) return 0;
  if (hasSymmetry(sym, rotatez2) && k%4) return 0;
  if (hasSymmetry(sym, rotatez3) && k%6) return 0;
  // if (hasSymmetry(sym, rotatey1 && ) return 0;
  // if (hasSymmetry(sym, reflectxz) && k==0 && 
  //    (j%4 && pi==0 || !(j%4) && pi==1)) return 0;
  
  return 1;
}


void extractSymmetryfromString(char** str, Symmetry* S)
{
  *S = 0;
  char *c, *s = *str;

  while (1) {
    c = strtok(s, ": >");
    if (c == NULL) {
      *str = c;
      return;
    }

    if (!strcmp(c, "parity0")) setSymmetry(S, parity0);
    else if (!strcmp(c, "parity1")) setSymmetry(S, parity1);
    else if (!strcmp(c, "spherical")) setSymmetry(S, spherical);
    else if (!strcmp(c, "axial0")) setSymmetry(S, axial0);
    else if (!strcmp(c, "axial1")) setSymmetry(S, axial1);
    else if (!strcmp(c, "axial2")) setSymmetry(S, axial2);
    else if (!strcmp(c, "axial3")) setSymmetry(S, axial3);
    else if (!strcmp(c, "axial4")) setSymmetry(S, axial4);
    else if (!strcmp(c, "rotatez2")) setSymmetry(S, rotatez2);
    else if (!strcmp(c, "rotatez3")) setSymmetry(S, rotatez3);
    else if (!strcmp(c, "rotatey1")) setSymmetry(S, rotatey1);
    else if (!strcmp(c, "reflectxz")) setSymmetry(S, reflectxz);
    else if (!strcmp(c, "reflectxy")) setSymmetry(S, reflectxy);
    else { 
      *str = c; 
      return;
    }
    s = NULL;
  }
}


// doesn't handle multiple symmetries yet
char* SymmetrytoStr(Symmetry S)
{
  char* str;

  str = malloc(40*sizeof(char));

  if (S==0) sprintf(str, "%s", "none");
  if (hasSymmetry(S, parity0))   sprintf(str, "%s", "parity0");
  if (hasSymmetry(S, parity1))   sprintf(str, "%s", "parity1");
  if (hasSymmetry(S, spherical)) sprintf(str, "%s", "spherical");
  if (hasSymmetry(S, axial))	 sprintf(str, "%s", "axial");
  if (hasSymmetry(S, axial0))    sprintf(str, "%s", "axial0");
  if (hasSymmetry(S, axial1))    sprintf(str, "%s", "axial1");
  if (hasSymmetry(S, axial2))    sprintf(str, "%s", "axial2");
  if (hasSymmetry(S, axial3))    sprintf(str, "%s", "axial3");
  if (hasSymmetry(S, axial4))    sprintf(str, "%s", "axial4");
  if (hasSymmetry(S, rotatez2))  sprintf(str, "%s", "rotatez2");
  if (hasSymmetry(S, rotatez3))  sprintf(str, "%s", "rotatez3");
  if (hasSymmetry(S, reflectxz)) sprintf(str, "%s", "reflectxy");
  if (hasSymmetry(S, reflectxy)) sprintf(str, "%s", "reflectxz");

  return str;
}
