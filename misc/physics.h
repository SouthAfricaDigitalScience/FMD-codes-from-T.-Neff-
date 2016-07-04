/**

  \file physics.h

  define Physical constants and properties


  (c) 2003 Thomas Neff

*/


#ifndef _PHYSICS_H
#define _PHYSICS_H


/// the names of the elements
extern const char* ElementName[]; 



extern double hbc;		///< hbarc [MeV fm]
extern double mproton;		///< proton mass [fm^-1]
extern double mneutron;		///< neutron mass [fm^-1]
extern double mnucleon;		///< nucleon mass [fm^-1]
extern double alpha;		///< finestructure constant	
extern double r2proton;		///< proton charge radius [fm^2]
extern double r2neutron;	///< neutron charge radius [fm^2]
extern double rho0;		///< nuclear matter density [fm^-3]
extern double muproton;		///< proton magnetic moment [mu_N]
extern double muneutron;	///< neutron magnetic moment [mu_N]
extern double gagv;             ///< beta-decay g_A/g_V


/// nucleon mass as a function of isospin
inline static double mass(int xi)
{
  if (xi == 1)
    return mproton;
  else if (xi == -1)
    return mneutron;

  return 0.0;
}


char* nucleusname(int A, int Z);

double r2charge(double r2p, int N, int Z);

double fproton(double q);
double fneutron(double q);

#endif
