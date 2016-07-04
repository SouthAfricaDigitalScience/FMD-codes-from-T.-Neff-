/**

  \file calctransitionsgt.c

  calculate GT transitions


  (c) 2004,2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/GamovTeller.h"
#include "fmd/Projection.h"

#include "numerics/zcw.h"
#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstatefin mcstateini"
	    "\n   -a             show all eigenstates\n", argv[0]);
    exit(-1);
  }

  int all=0;
  int hermit=0;

  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "a")) != -1)
    switch (c) {
    case 'a':
      all=1;
      break;
    }

  char* mcstatefilefin = argv[optind];
  char* mcstatefileini = argv[optind+1];
  char** mbfilefin; char** mbfileini;

  // open multiconfigfile
  Projection P;
  SlaterDet *Qfin, *Qini;
  Symmetry *Sfin, *Sini;
  Eigenstates Efin, Eini;
  int nfin, nini;

  readMulticonfigfile(mcstatefileini, &mbfileini, &P, &Qini, &Sini, &Eini, &nini);
  readMulticonfigfile(mcstatefilefin, &mbfilefin, &P, &Qfin, &Sfin, &Efin, &nfin);

  int direction;
  if (Qfin[0].Z - Qini[0].Z == 1) {
    fprintf(stderr, "GT+ transitions");
    direction = +1;
  }
  if (Qfin[0].Z - Qini[0].Z == -1) {
    fprintf(stderr, "GT- transitions");
    direction = -1;
  }

 
  ManyBodyOperator OpGT;

  if (direction == +1)
    OpGT = OpGTplus;
  if (direction == -1)
    OpGT = OpGTminus;


  int a,b; 

  void** gtme[nfin*nini];

  for (b=0; b<nini; b++)	
    for (a=0; a<nfin; a++)
      gtme[a+b*nfin] = initprojectedMBME(&P, &OpGT);


  // read or calculate matrix elements
  for (b=0; b<nini; b++)
    for (a=0; a<nfin; a++)
      if (readprojectedMBMEfromFile(mbfilefin[a], mbfileini[b], &P,
                                    &OpGT, Sfin[a], Sini[b], gtme[a+b*nfin])) {
	calcprojectedMBME(&P, &OpGT, &Qfin[a], &Qini[b], Sfin[a], Sini[b],
                          gtme[a+b*nfin]);
        writeprojectedMBMEtoFile(mbfilefin[a], mbfileini[b], &P,
                                 &OpGT, Sfin[a], Sini[b], gtme[a+b*nfin]);
      } 
  
  fprintf(stderr, "calculate transition strengths\n");

  void**** gttrans = initprojectedtransitionVector(&P, &OpGT, nfin, nini);
  calctransitionprojectedMBME(&P, &OpGT, gtme, Sfin, Sini, &Efin, &Eini, gttrans);

  // output

  char outfile[511];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s--%s.transgt", 
           stripstr(mcstatefilefin, tostrip), stripstr(mcstatefileini, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  if (direction == +1)
    writeprojectedtransitionGTplus(outfp, &P, gttrans, &Efin, &Eini);
  else
    writeprojectedtransitionGTminus(outfp, &P, gttrans, &Efin, &Eini);

  fclose(outfp);

 cleanup:

  return 0;
}


