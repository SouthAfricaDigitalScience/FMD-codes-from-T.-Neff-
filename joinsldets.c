/**

  \file joinnuclei.c

  read two parameterizations and create joined parameter file

  
  (c) 2004 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "fmd/Parameterization.h"
#include "fmd/ParameterizationFMD.h"
#include "fmd/SlaterDet.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  int cm=1;
  double d=10.0;
  double E=0.0, p=0.0;
  double b=0.0;

  double xA[3] = {0.0, 0.0, 0.0}, xB[3] = {0.0, 0.0, 0.0};
  double vA[3] = {0.0, 0.0, 0.0}, vB[3] = {0.0, 0.0, 0.0};
  double alphaA=0.0, alphaB=0.0;
  double betaA=0.0, betaB=0.0;
  double gammaA=0.0, gammaB=0.0;

  int c, i;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s [OPTIONS] NucA NucB JoinedNuc\n"
            "\n  -c            center of mass system (default)"
            "\n  -l            lab system, Nucleus A at x=0 in rest"
            "\n  -d dist       distance of nuclei [fm]"
	    "\n  -p momentum   relative momentum of nuclei [fm^-1]"
            "\n  -b impact     impact parameter (default 0 fm)"
            "\n  -r Aalpha,beta,gamma     rotate NucA by alpha, beta, gamma (degrees)"
            "\n  -r Balpha,beta,gamma     rotate NucB by alpha, beta, gamma (degrees)"
	    "\n  -E energy     E/A in lab frame (default 0 MeV)\n",
    argv[0]);
    exit(-1);
  }

  /* manage command line options */

  while ((c = getopt(argc, argv, "cld:p:E:b:r:")) != -1)
    switch (c) {
    case 'c':
      cm=1;
      break;
    case 'l':
      cm=0;
      break;
    case 'd':
      d = atof(optarg);
      break;
    case 'p':
      p = atof(optarg);
      break;
    case 'b':
      b = atof(optarg);
      break;
    case 'E':
      E = atof(optarg);
      break;
    case 'r':
      if (optarg[0] == 'A') {
	sscanf(++optarg, "%lf,%lf,%lf", &alphaA, &betaA, &gammaA);
	alphaA *= M_PI/180;
	betaA *= M_PI/180;
	gammaA *= M_PI/180;
      }
      else if (optarg[0] == 'B') {
	sscanf(++optarg, "%lf,%lf,%lf", &alphaB, &betaB, &gammaB);
	alphaB *= M_PI/180;
	betaB *= M_PI/180;
	gammaB *= M_PI/180;
      } 
      break;
    }
     
 
  char* QAfile = argv[optind];
  char* QBfile = argv[optind+1];
  char* Qfile = argv[optind+2];

  SlaterDet QA, QB, Q;
  if (readSlaterDetfromFile(&QA, QAfile)) {
    fprintf(stderr, "couldn't read %s\n", QAfile);
    exit(-1);
  }
  if (readSlaterDetfromFile(&QB, QBfile)) {
    fprintf(stderr, "couldn't read %s\n", QBfile);
    exit(-1);
  }

  double massA, massB, mass;

  massA = QA.Z*mproton+QA.N*mneutron;
  massB = QB.Z*mproton+QB.N*mneutron;
  mass = massA+massB;

  if (E != 0.0) {
    E = QA.A*E/hbc;
    p = sqrt(2.0*massA*E);
  }

  if (cm) {
    xA[2] = -massB/mass*d;
    xB[2] = massA/mass*d;
    xA[1] = -massB/mass*b;
    xB[1] = massA/mass*b;

    // vA[2] = p/massA*(1-massA/mass);
    // vB[2] = -p/massA*(massA/mass);
    vA[2] = -p/massA;
    vB[2] = p/massB;
  } else {
    // xA[2] = -d; xA[1] = -b;
    // vA[2] = p/massA;
    xB[2] = d; xB[1] = b;
    vB[2] = p/massB;
  }

  rotateSlaterDet(&QA, alphaA, betaA, gammaA);
  boostSlaterDet(&QA, vA);
  moveSlaterDet(&QA, xA);

  rotateSlaterDet(&QB, alphaB, betaB, gammaB);
  boostSlaterDet(&QB, vB);
  moveSlaterDet(&QB, xB);

  joinSlaterDets(&QA, &QB, &Q);
  
  Parameterization P = ParameterizationFMD;
  Para q;
  SlaterDetinitFMD(&Q, &q);

  // write joined SlaterDet parameters
  
  FILE* outfp;
 
  fprintf(stderr, "... writing joined SlaterDet to file %s\n", Qfile);
  if (!(outfp = fopen(Qfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", Qfile);
    exit(-1);
  }

  fprintinfo(outfp);

  fprintf(outfp, "\n# joined SlaterDets from %s and %s\n",
	  QAfile, QBfile);
  
  fprintf(outfp, "\n# Parameterization\n");
  fprintf(outfp, "<Parameterization %s>\n", P.name);
  P.Parawrite(outfp, &q);

  fprintf(outfp, "\n# SlaterDet\n");
  writeSlaterDet(outfp, &Q);

  fclose(outfp);

  exit(0);
}
