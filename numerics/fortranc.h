/*
  fortranc.h

  naming conventions for fortran symbols

  (c) Thomas Neff 2003

*/

#ifndef _FORTRANC_H
#define _FORTRANC_H


#ifdef __aix__
#define FORTRAN(name) name
#endif

#ifdef __linux__
#define FORTRAN(name) name##_
#endif


#endif
