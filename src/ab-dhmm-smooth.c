/*    File:  ab-dhmm-smooth.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Smooth the results of an approximate Bayes chain.
 *           (Should also be usable for a penalized EM-step chain.)
 *
 *   Usage:  ab-dhmm-smooth posterior data smooth
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
  int I, M, T;     /* states, observation types, sequence length */
  int *X;          /* integer array */
  int i, t;        /* indices */

  double *scale, **alpha, **beta, **gamma, ***xi;

  DHMM prior, star;
  FILE *DATA, *PRIOR, *SMOOTH;

  /* Initialize */
  /* To do:  1. Put in argument-checking code here. */
  
  /* read in the prior */
  PRIOR = fopen(argv[1], "r");
  read_DHMM(PRIOR, &prior);
  fclose(PRIOR);  

  /* read in the data */
  DATA = fopen(argv[2], "r");
  read_sequence(DATA, &T, &X);
  fclose(DATA);

  /* Debug code */
  /* fprintf(stdout, "T: %d and X[T]: %d\n", T, X[T]); */

  /* allocate working space */
  I = prior.I;
  M = prior.M;

  alpha = array_double_2d(1, I, 1, T);
  beta  = array_double_2d(1, I, 1, T);
  gamma = array_double_2d(1, I, 1, T);
  xi    = array_double_3d(1, I, 1, I, 1, T);

  new_DHMM(&star, I, M);                   /* initialize the star matrix */
  scale = array_double_1d(1, T);           /* initialize scaling array */

  /* Process */
  
  /* We use the fact that the normalized gamma and xi arrays
   * contain respectively the posterior probabilities of
   * states and observations.
   */

  /* M-step */
  Calculate_star(&prior, &star);

  /* E-step */
  /* We do not want to rescale the xi and gamma arrays here */
  Forward_Scaled(&star, T, X, alpha, scale);
  Backward_Scaled(&star, T, X, beta, scale);

  Accumulate_gamma(&star, T, X, alpha, beta, gamma);

  /* To do:  1.  Put the above and the below into separate utility
   *             routines.
   */

  /* write results */
  SMOOTH = fopen(argv[3], "w");

  for (t=1; t<=T; t++) {
    for (i=1; i<=I; i++) {
      fprintf(SMOOTH, "%12.3f ", gamma[i][t]);
    }
    fprintf(SMOOTH, "\n");
  }

  fclose(SMOOTH);

  /* Clean-up */

  free_array_double_1d(scale, 1, T);
  free_array_double_2d(alpha, 1, I, 1, T);
  free_array_double_2d( beta, 1, I, 1, T);
  free_array_double_2d(gamma, 1, I, 1, T);
  free_array_double_3d(   xi, 1, I, 1, I, 1, T);
  free_array_int_1d(    X, 1, T);
  free_DHMM(&prior);
  free_DHMM(&star);
}
