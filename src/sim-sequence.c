/*    File:  sim-sequence.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Simulate a realization from a DHMM.
 *
 *    Note:  1.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   Usage:  sim-sequence dhmm T
 */

static char rcsid[] = "$Id: programs.nw,v 1.10 1999/07/26 08:51:31 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "dhmm.h"               /* uses FILE so we need <stdio.h> first */
#include "utilities.h"

/*********************************************************************
 * PROGRAM 
 ********************************************************************/

void main (int argc, char **argv, char **env)
{
  /* Initialize */
  int T;
  int *X;

  DHMM prior;
  FILE *PRIOR;

  unsigned long seed;      /* pseudorandom generator seed */

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */
  
  set_seed(&seed, 0);

  /* read in the DHMM */
  PRIOR = fopen(argv[1], "r");
  read_DHMM(PRIOR, &prior);
  fclose(PRIOR);  

  /* get the sequence length */
  T = atoi(argv[2]);

  /* generate a sequence from the DHMM realization */
  generate_sequence(&seed, &prior, T, &X);

  /* write out the sequence */
  write_sequence(stdout, T, X);

  /* Clean-up */
  save_seed(&seed);

  free_DHMM(&prior);
}





