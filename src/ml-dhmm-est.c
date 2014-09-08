/*    File:  ml-dhmm-est.c
 *
 *  Author:  Jerome V. Braun <jerome.braun@kmri.com>
 * Company:  Kings Mountain Research
 *         
 * Purpose:  Carry out the maximum likelihood estimation
 *           of parameters of the hidden Markov model using
 *           the Baum-Welch (EM) algorithm.
 *
 *    Note:  Only handles one sequence.
 *
 *   Usage:  ml-dhmm-est initial data estimate
 */

static char rcsid[] = "$Id: programs.nw,v 1.6 1999/06/14 03:47:43 Jerome Braun Exp $";

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
  int I, M, T;

  double **alpha, **beta, **gamma, ***xi;
  int *X;

  DHMM  dhmm;
  FILE *INITIAL, *DATA, *ESTIMATE;

  unsigned long seed;            /* pseudorandom generator seed */
  
  /* Set-up */
  /* To do:  1. Put in argument-checking code here. */

  /* Initialize */
  
  INITIAL = fopen(argv[1], "r");   /* read in the initial estimate */
    read_DHMM(INITIAL, &dhmm);
  fclose(INITIAL);  

  DATA = fopen(argv[2], "r");    /* read in the data */
    read_sequence(DATA, &T, &X);
  fclose(DATA);

  /* allocate working space */
  I = dhmm.I;
  M = dhmm.M;

  alpha = array_double_2d(1, I, 1, T);
  beta  = array_double_2d(1, I, 1, T);
  gamma = array_double_2d(1, I, 1, T);
  xi    = array_double_3d(1, I, 1, I, 1, T);

  /* Process */

  /* run the maximum likelihood algorithm and write results */

  Baum_Welch(&dhmm, T, X, alpha, beta, gamma, xi);

  ESTIMATE = fopen(argv[3], "w");
  write_DHMM(ESTIMATE, &dhmm);
  fclose(ESTIMATE);

  /* Clean-up */

  free_array_double_2d(alpha, 1, I, 1, T);
  free_array_double_2d( beta, 1, I, 1, T);
  free_array_double_2d(gamma, 1, I, 1, T);
  free_array_double_3d(   xi, 1, I, 1, I, 1, T);
  free_array_int_1d(    X, 1, T);
  free_DHMM(&dhmm);
}

