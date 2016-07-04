/**

  \file physics.c

  define Physical constants and properties


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "physics.h"


const char* ElementName[] = 
  {"n", "H", "He",
   "Li", "Be", "B", "C", "N", "O", "F", "Ne",
   "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
   "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
   "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
   "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
   "Pd", "Ag", "Cd", "In"};


// Standard values taken from PDB

double hbc = 197.327053;
double mproton = 938.27231/197.327053;	
double mneutron = 939.56563/197.327053;	
double mnucleon = 0.5*(938.27231+939.56563)/197.327053;	
double alpha = 1/137.0359895;		

double r2proton = 0.801025;     // Rp = 0.895(18) fm, I. Sick, PLB 576 (2003) 62
double r2neutron = -0.120; 	// Rn^2 = -0.120(5) fm^2, S. Kopecky, PRC 56 (1997) 2229

double rho0 = 0.17;

double muproton = 2.792847337;
double muneutron = -1.9130427;

double gagv = -1.2695;           // -1.2695 \pm 0.0029 


char* nucleusname(int A, int Z)
{
  char* name = malloc(5*sizeof(char));
  sprintf(name, "%s%d", ElementName[Z], A);
  
  return name;
}


double r2charge(double r2p, int N, int Z)
{
  return (r2p + r2proton + 1.0*N/Z*r2neutron);
}


// proton and neutron formfactor taken from 
// Z. Phys. A - Atomic Nuclei 323, 13-25 (1986)

/*
double fproton(double q)
{
  const double a[4] = { 0.312, 1.312, -0.709, 0.085};
  const double b[4] = { 0.16667, 0.06658, 0.02269, 0.006485};
  double f = 0.0;
  int i;

  for (i=0; i<4; i++)
    f += a[i]/(1+b[i]*q*q);

  return f;
}


double fneutron(double q)
{
  const double a[2] = { 1, -1};
  const double b[2] = { 0.04833, 0.05833};
  double f = 0.0;
  int i;

  for (i=0; i<2; i++)
    f += a[i]/(1+b[i]*q*q);

  return f;
}
*/

// proton and neutron formfactor taken from
// Friedrich, Walcher, Eur. Phys. J. A 17 607-623 (2003)
// phenomenological Fit, Table 2

#define SQR(x) ((x)*(x))

double fproton(double q)
{
  double a10 = 1.041, a11 = 0.765*SQR(1000/hbc);
  double a20 = -0.041, a21 = 6.2*SQR(1000/hbc);
  double ab = -0.23/SQR(1000/hbc), Qb = 0.07*(1000/hbc), 
    sigb = 0.27*(1000/hbc);

  return (a10/SQR(1+SQR(q)/a11) + a20/(SQR(1+SQR(q)/a21)) +
	  ab*SQR(q)*(exp(-0.5*SQR((q-Qb)/sigb))+exp(-0.5*SQR((q+Qb)/sigb))));
}

double fneutron(double q)
{
  double a10 = 1.04, a11 = 1.73*SQR(1000/hbc);
  double a20 = -1.04, a21 = 1.54*SQR(1000/hbc);
  double ab = 0.23/SQR(1000/hbc), Qb = 0.29*(1000/hbc), 
    sigb = 0.20*(1000/hbc);

  return (a10/SQR(1+SQR(q)/a11) + a20/(SQR(1+SQR(q)/a21)) +
	  ab*SQR(q)*(exp(-0.5*SQR((q-Qb)/sigb))+exp(-0.5*SQR((q+Qb)/sigb))));
}



