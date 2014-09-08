/*    File:  pml-dhmm-sim.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Simulate a realization from a model drawn from 
 *           the posterior mode model.
 *
 *   To do:  1. Consider sending the output to a file.
 *
 *    Note:  1.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   Usage:  pml-dhmm-est posterior_mode T
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
  int T;            /* sequence length */
  int *X;           /* integer array */

  DHMM posterior_mode;
  FILE *POSTERIOR_MODE;

  unsigned long seed;            /* pseudorandom generator seed */

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */
  
  set_seed(&seed, 1);      /* default 4 * 0 + 1 */
  
  /* get the sequence length to simulate */
  T = atoi(argv[2]);  

  /* read in the posterior mode */
  POSTERIOR_MODE = fopen(argv[2], "r");
  read_DHMM(POSTERIOR_MODE, &posterior_mode);
  fclose(POSTERIOR_MODE);  

  /* generate a sequence from the posterior mode */
  generate_sequence(&seed, &posterior_mode, T, &X);

  /* write out the sequence */
  write_sequence(stdout, T, X);

  /* Clean-up */
  save_seed(&seed);

  free_array_int_1d(    X, 1, T);
  free_DHMM(&posterior_mode);
}

