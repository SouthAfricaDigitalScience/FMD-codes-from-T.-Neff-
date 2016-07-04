/**

  \file gnd2fmdpara.c

  convert old gnd format into FMD parameterization


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "misc/utils.h"
#include "misc/physics.h"

int main(int argc, char *argv[])
{
  createinfo(argc, argv);

  int c;
  int spinfree=0;
  int multi=0;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s gndfile fmdparafile"
	    "\n		-s 	free spin gndfile"
	    "\n		-m	multi gauss gndfile\n",
	    argv[0]);
    exit(-1);
  }
  
  while((c = getopt(argc, argv, "sm")) != -1)
    switch (c) {
    case 's':
      spinfree = 1;
      break;
    case 'm':
      multi = 1;
      break;
    }
  
  char* gndfile = argv[optind];
  char* fmdparafile = argv[optind+1];

  char name[255];
  int A; int ng;
  char buf[1024];
  
  FILE* gndfp = fopen(gndfile, "r");

  do
    fgets(buf, 255, gndfp);
  while (buf[0] == '#' || buf[0] == '\n');
  sscanf(buf, "%s", name);

  do
    fgets(buf, 255, gndfp);
  while (buf[0] == '#' || buf[0] == '\n');
  if (multi)
    sscanf(buf, "%d %d", &A, &ng);
  else {
    sscanf(buf, "%d", &A);
    ng = 1;
  }

  int i;

  do
    fgets(buf, 255, gndfp);
  while (buf[0] == '#' || buf[0] == '\n');

  int xi[A*ng]; int chi; 
  double alpha, beta;
  double cre, cim;
  double are, aim, x[3], p[3];
  complex double newA[A*ng], newB[A*ng][3];
  complex double newChi[A*ng][2];
  int j;

  for (i=0; i<A*ng; i++) {
    if (spinfree) {
      sscanf(buf, "%d  %lf %lf %lf  %lf %lf %lf  %lf %lf %lf %lf", 
	     &xi[i], &x[0], &x[1], &x[2], &p[0], &p[1], &p[2], &are, &aim, 
	     &alpha, &beta);

      newChi[i][0] = cos(0.5*alpha);
      newChi[i][1] = sin(0.5*alpha)*cexp(I*beta);
    } else if (multi) {
      sscanf(buf, "%lf %lf  %d %d  %lf %lf %lf  %lf %lf %lf  %lf %lf",
	     &cre, &cim,
	     &xi[i], &chi, &x[0], &x[1], &x[2], 
	     &p[0], &p[1], &p[2], &are, &aim);
      if (chi == +1) { newChi[i][0] = cre+I*cim; newChi[i][1] = 0.0; }
      else           { newChi[i][0] = 0.0; newChi[i][1] = cre+I*cim; }
    } else {
      sscanf(buf, "%d  %d %lf %lf %lf  %lf %lf %lf  %lf %lf", 
	     &xi[i], &chi, &x[0], &x[1], &x[2], 
	     &p[0], &p[1], &p[2], &are, &aim);
      if (chi == +1) { newChi[i][0] = 1.0; newChi[i][1] = 0.0; }
      else           { newChi[i][0] = 0.0; newChi[i][1] = 1.0; }
    }
    newA[i] = are + I*aim;
    for (j=0; j<3; j++)
      newB[i][j] = x[j] + I*newA[i]* p[j];
    
    fgets(buf, 255, gndfp);
  } 
  fclose(gndfp);

  int Z=0, N=0;
  for (i=0; i<A; i++) {
    if (xi[ng*i] == +1) Z++;
    if (xi[ng*i] == -1) N++;
  }

  // write FMD parameters to file

  FILE* parafp = fopen(fmdparafile, "w");

  fprintinfo(parafp);

  fprintf(parafp, "\n<Parameterization FMD>\n");
  fprintf(parafp, "<Para FMD %s>\n", name);
  fprintf(parafp, "%2d %2d %2d\n", A, Z, N);
  fprintf(parafp, "%2d\n", A*ng);
  for (i=0; i<A*ng; i++)
    fprintf(parafp, "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	    i / ng, xi[i], 
	    creal(newChi[i][0]), cimag(newChi[i][0]), 
	    creal(newChi[i][1]), cimag(newChi[i][1]),
	    creal(newA[i]),      cimag(newA[i]), 
	    creal(newB[i][0]),   cimag(newB[i][0]),
	    creal(newB[i][1]),   cimag(newB[i][1]), 
	    creal(newB[i][2]),   cimag(newB[i][2]));
  fprintf(parafp, "</Para>\n");

  fprintf(parafp, "\n\n<SlaterDet>\n");
  fprintf(parafp, "%2d %2d %2d\n", A, Z, N);
  fprintf(parafp, "%2d\n", A*ng);
  for (i=0; i<A*ng; i++)
    fprintf(parafp, "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	    i / ng, xi[i], 
	    creal(newChi[i][0]), cimag(newChi[i][0]), 
	    creal(newChi[i][1]), cimag(newChi[i][1]),
	    creal(newA[i]),      cimag(newA[i]), 
	    creal(newB[i][0]),   cimag(newB[i][0]),
	    creal(newB[i][1]),   cimag(newB[i][1]), 
	    creal(newB[i][2]),   cimag(newB[i][2]));
  fprintf(parafp, "</SlaterDet>\n");

  fclose(parafp);

  return 0;
}
