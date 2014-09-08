/*    File:  ab-dhmm-est.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Carry out the approximate Bayes estimation
 *           of the posterior distribution.
 *
 *    Note:  1.  Only handles one sequence.
 *           2.  Relies on a .SEED file to initialize the pseudorandom
 *               number generator seed.  If there is none, the start value
 *               4 * 0 + 1 = 1 is used and a .SEED file is created.
 *
 *   Usage:  ab-dhmm-est prior data posterior
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

  int I, M, T;             /* states, observations, sequence length */
  int *X;                  /* integer array */

  double **alpha, **beta, **gamma, ***xi;

  DHMM posterior, prior;
  FILE *PRIOR, *DATA, *POSTERIOR;

  unsigned long seed;      /* pseudorandom generator seed */
  
  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */

  set_seed(&seed, 0);
  
  /* read in the prior */
  PRIOR = fopen(argv[1], "r");
  read_DHMM(PRIOR, &prior);
  fclose(PRIOR);  

  /* read in the data */
  DATA = fopen(argv[2], "r");
  read_sequence(DATA, &T, &X);
  fclose(DATA);

  /* allocate working space */
  I = prior.I;
  M = prior.M;

  new_DHMM(&posterior, I, M);

  alpha = array_double_2d(1, I, 1, T);
  beta  = array_double_2d(1, I, 1, T);
  gamma = array_double_2d(1, I, 1, T);
  xi    = array_double_3d(1, I, 1, I, 1, T);

  /* Process */

  /* we must initialize the posterior mode or the algorithm */
  /* will not work correctly */
  initialize_DHMM(&seed, &posterior);

  /* run the approximate Bayes algorithm and write results */
  Approximate_Bayes(&posterior, &prior, T, X, alpha, beta, gamma, xi);

  /* write posterior mode to a file */
  POSTERIOR = fopen(argv[3], "w");
  write_DHMM(POSTERIOR, &posterior);
  fclose(POSTERIOR);

  /* Clean-up */
  save_seed(&seed);

  free_array_double_2d(alpha, 1, I, 1, T);
  free_array_double_2d( beta, 1, I, 1, T);
  free_array_double_2d(gamma, 1, I, 1, T);
  free_array_double_3d(   xi, 1, I, 1, I, 1, T);
  free_array_int_1d(    X, 1, T);
  free_DHMM(&posterior);
  free_DHMM(&prior);
}
