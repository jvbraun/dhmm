/*    File:  initialize-dhmm.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Build, randomly initialize, and write out a DHMM structure.
 *
 *    Note:  1.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   To do:  1. Consider sending the output to a file.
 *
 *   Usage:  initialize-dhmm I M
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
  /* Start-up */

  int I, M;             /* states, observations */

  DHMM dhmm;

  unsigned long seed;   /* pseudorandom generator seed */

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */

  set_seed(&seed, 0);

  I = atoi(argv[1]);
  M = atoi(argv[2]);

  new_DHMM(&dhmm, I, M);

  /* initialize the DHMM realization */
  initialize_DHMM(&seed, &dhmm);

  /* write out the sequence */
  write_DHMM(stdout, &dhmm);

  /* Clean-up */
  save_seed(&seed);

  free_DHMM(&dhmm);
}
