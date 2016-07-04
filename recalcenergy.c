/**

  \file recalcenergy.c

  read parameterization and
  calculate energy of Slater determinant


  (c) 2003 Thomas Neff

*/

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/Parameterization.h"
#include "fmd/ParameterizationFMD.h"
#include "fmd/SlaterDet.h"
#include "fmd/Interaction.h"
#include "fmd/Hamiltonian.h"
#include "fmd/Observables.h"

#include "misc/physics.h"
#include "misc/utils.h"


int main(int argc, char *argv[])
{
  createinfo(argc, argv);

  int cm=1;
  int overwrite=0;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction parafile\n"
	    "\n   -o              overwrite parameter file\n",
	    argv[0]);
    exit(-1);
  }
  
  int c;
  while((c = getopt(argc, argv, "o")) != -1)
    switch (c) {
    case 'o':
      overwrite = 1;
      break;
    }
  
  char* interactionfile = argv[optind];
  char* parafile = argv[optind+1];

  Interaction Int;
  readInteractionfromFile(&Int, interactionfile);
  Int.cm = cm;

  Parameterization P;
  Para q;
  readParafromFile(&P, &q, parafile); 

  SlaterDet Q;
  SlaterDetAux X;

  P.ParainitSlaterDet(&q, &Q);
  initSlaterDetAux(&Q, &X);
  P.ParatoSlaterDet(&q, &Q);

  calcSlaterDetAux(&Q, &X);

  // moveboostSlaterDet(&Q, &X);
  moveboostorientSlaterDet(&Q, &X);

  // normalize SlaterDet
  normalizeSlaterDet(&Q, &X);

  // calculate the observables
  Observables Obs;
  
  calcSlaterDetAux(&Q, &X);
  calcObservables(&Int, &Q, &Obs);

  FILE* outfp;
  if (overwrite) {

    backup(parafile);
    fprintf(stderr, "... writing Parameters to file %s\n", parafile);
    if (!(outfp = fopen(parafile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", parafile);
      exit(-1);
    }
  } else {
    outfp = stdout;
  }

  fprintinfo(outfp);

  fprintf(outfp, "\n# calculated %s for %s in %s parameterization\n"
	  "# using %s interaction\n", 
	  cm ? "< Hintr >" : "< H >", q.name, P.name, Int.name);

  fprintObservables(outfp, &Int, &Q, &Obs);

  fprintf(outfp, "\n# Parameterization\n");
  fprintf(outfp, "<Parameterization %s>\n", P.name);
  P.Parawrite(outfp, &q);
 
  fprintf(outfp, "\n# SlaterDet\n");
  writeSlaterDet(outfp, &Q);

  if (overwrite) 
    fclose(outfp);

  return 0;
}
