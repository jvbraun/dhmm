/*    File:  sim-dhmm.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Simulate a realization from a prior on DHMM..
 *
 *    Note:  1.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   Usage:  sim-dhmm prior
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

  DHMM prior, dhmm;
  FILE *PRIOR;

  unsigned long seed;      /* pseudorandom generator seed */

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */
  
  set_seed(&seed, 0);

  /* read in the prior */
  PRIOR = fopen(argv[1], "r");
  read_DHMM(PRIOR, &prior);
  fclose(PRIOR);  

  /* allocate working space for the DHMM realization */
  new_DHMM(&dhmm, prior.I, prior.M);

  /* generate the DHMM realization */
  DHMM_deviate(&seed, &prior, &dhmm);

  /* write out the DHMM */
  write_DHMM(stdout, &dhmm);

  /* Clean-up */
  save_seed(&seed);

  free_DHMM(&dhmm);
  free_DHMM(&prior);
}

