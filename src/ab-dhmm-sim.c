/*    File:  ab-dhmm-sim.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Simulate a realization from a model drawn from 
 *           the distribution on chains.
 *
 *    Note:  1.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   To do:  1. Consider sending the output to a file.
 *
 *   Usage:  ab-dhmm-est prior T
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

  int I, M, T;             /* states, observations, sequence length */
  int *X;                  /* integer array */

  DHMM prior, dhmm;
  FILE *PRIOR;

  unsigned long seed;      /* pseudorandom generator seed */

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */
  
  set_seed(&seed, 0);

  /* get the sequence length to simulate */
  T = atoi(argv[2]);  

  /* read in the prior */
  PRIOR = fopen(argv[1], "r");
    read_DHMM(PRIOR, &prior);
  fclose(PRIOR);  

  /* allocate working space for the DHMM realization */
  I = prior.I;
  M = prior.M;

  new_DHMM(&dhmm, I, M);

  /* generate the DHMM realization */
  DHMM_deviate(&seed, &prior, &dhmm);

  /* generate a sequence from the DHMM realization */
  generate_sequence(&seed, &dhmm, T, &X);

  /* write out the sequence */
  write_sequence(stdout, T, X);

  /* Clean-up */
  save_seed(&seed);

  free_array_int_1d(    X, 1, T);
  free_DHMM(&dhmm);
  free_DHMM(&prior);
}

