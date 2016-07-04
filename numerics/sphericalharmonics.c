/**

   \file sphericalharmonics.c

   gives Spherical Harmonics Y(l,m,theta,phi)
   brute force no brain table based implementation

   l is 0,2,4,6,8
   
   (c) 2004 Thomas Neff
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "sphericalharmonics.h"


#define SQR(x)	((x)*(x))
#define CBE(x)	((x)*(x)*(x))
#define QRT(x)  ((x)*(x)*(x)*(x))

complex double Y(int l, int m, double theta, double phi)
{
  switch (l) {
  case 0 :
    return (0.5/sqrt(M_PI));
  case 2 :
    switch (m) {
    case +2 :
      return (-0.5*cexp(I*phi)*sqrt(1.5/M_PI)*sin(theta));
    case 0 :
      return (0.5*sqrt(3.0/M_PI)*cos(theta));
    case -2 :
      return (0.5*cexp(-I*phi)*sqrt(1.5/M_PI)*sin(theta));
    }
  case 4 :
    switch (m) {
    case +4 :
      return (0.25*cexp(2*I*phi)*sqrt(7.5/M_PI)*SQR(sin(theta)));
    case +2 :
      return (-0.5*cexp(I*phi)*sqrt(7.5/M_PI)*cos(theta)*sin(theta));
    case 0 :
      return (0.25*sqrt(5/M_PI)*(3*SQR(cos(theta))-1));
    case -2 :
      return (0.5*cexp(-I*phi)*sqrt(7.5/M_PI)*cos(theta)*sin(theta));
    case -4 :
      return (0.25*cexp(-2*I*phi)*sqrt(7.5/M_PI)*SQR(sin(theta)));
    }
  case 6 :
    switch (m) {
    case +6 :
      return (-0.125*cexp(3*I*phi)*sqrt(35/M_PI)*CBE(sin(theta)));
    case +4 :
      return (0.25*cexp(2*I*phi)*sqrt(52.5/M_PI)*cos(theta)*SQR(sin(theta)));
    case +2 :
      return (-0.125*cexp(I*phi)*sqrt(21/M_PI)*(5*SQR(cos(theta))-1)*sin(theta));
    case 0 :
      return (0.25*sqrt(7/M_PI)*(5*CBE(cos(theta))-3*cos(theta)));
    case -2 :
      return (0.125*cexp(-I*phi)*sqrt(21/M_PI)*(5*SQR(cos(theta))-1)*sin(theta));
    case -4 :
      return (0.25*cexp(-2*I*phi)*sqrt(52.5/M_PI)*cos(theta)*SQR(sin(theta)));
    case -6 :
      return (0.125*cexp(-3*I*phi)*sqrt(35/M_PI)*CBE(sin(theta)));
    }
  case 8 :
    switch (m) {
    case +8 :
      return (0.1875*cexp(4*I*phi)*sqrt(17.5/M_PI)*QRT(sin(theta)));
    case +6 :
      return (-0.375*cexp(3*I*phi)*sqrt(35/M_PI)*cos(theta)*CBE(sin(theta)));
    case +4 :
      return (0.375*cexp(2*I*phi)*sqrt(2.5/M_PI)*(7*SQR(cos(theta))-1)*SQR(sin(theta)));
    case +2 :
      return (-0.375*cexp(I*phi)*sqrt(5/M_PI)*(7*SQR(cos(theta))-3)*cos(theta)*sin(theta));
    case 0 :
      return (0.1875/sqrt(M_PI)*(35*QRT(cos(theta))-30*SQR(cos(theta))+3));
    case -2 :
      return (0.375*cexp(-I*phi)*sqrt(5/M_PI)*(7*SQR(cos(theta))-3)*cos(theta)*sin(theta));
    case -4 :
      return (0.375*cexp(-2*I*phi)*sqrt(2.5/M_PI)*(7*SQR(cos(theta))-1)*SQR(sin(theta)));
    case -6 :
      return (0.375*cexp(-3*I*phi)*sqrt(35/M_PI)*cos(theta)*CBE(sin(theta)));
    case -8 :
      return (0.1875*cexp(-4*I*phi)*sqrt(17.5/M_PI)*QRT(sin(theta)));
    }
  default:
    fprintf(stderr, "sphericalharmonics.c: Y(l,m,theta,phi) not implemented for l>8\n");
    exit(-1);
  }
}
