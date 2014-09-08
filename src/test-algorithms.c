/*      Program:  test-algorithms.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Test the discrete hidden Markov model algorithms.
 */

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

void main (void)
{
  int I, M, T;

  double **alpha, **beta, **gamma, ***xi;
  int *X;

  DHMM dhmm, model, prior;
  FILE *fp;
  unsigned long seed;  /* pseudorandom generator seed */
  
  /* Set up the environment */

  set_seed(&seed, 4); /* set the pseudorandom generator seed */ 

  I = 2;    /* number of states */
  M = 2;    /* number of observation types */
  T = 100;  /* sequence length */

  /* allocate the DHMM's and working space */
  new_DHMM(&dhmm, I, M);    /* estimate */
  new_DHMM(&model, I, M);   /* true model */
  new_DHMM(&prior, I, M);   /* prior */

  alpha = array_double_2d(1, I, 1, T);
  beta  = array_double_2d(1, I, 1, T);
  gamma = array_double_2d(1, I, 1, T);
  xi    = array_double_3d(1, I, 1, I, 1, T);

  /* set the prior */
  prior.A[1][1] = 1.0; prior.A[1][2] = 2.0;
  prior.A[2][1] = 4.0; prior.A[2][2] = 1.0;

  prior.B[1][1] = 1.0; prior.B[1][2] = 3.0;
  prior.B[2][1] = 2.0; prior.B[2][2] = 1.0;

  prior.pi[1] = 1.0;   prior.pi[2] = 2.0;

  /* generate a hidden Markov model from the prior */
  DHMM_deviate(&seed, &prior, &model);

  /* generate a sequence from the hidden Markov model */
  generate_sequence(&seed, &model, T, &X);   /* allocates space for X */

  /* initialize the hidden Markov estimate randomly */
  initialize_DHMM(&seed, &dhmm);

  /* Display the information so far */
  fp = fopen("t-a.out", "w");
  fprintf(fp, "* This output shows the result of running the different discrete\n");
  fprintf(fp, "* hidden Markov model estimation algorithms.\n\n");
  fprintf(fp, "* The prior distribution for the hidden Markov model:\n");
  write_DHMM(fp, &prior);
  fprintf(fp, "\n* A realization from the prior distribution:\n");
  write_DHMM(fp, &model);
  fprintf(fp, "\n* A sequence generated from the realization:\n");
  write_sequence(fp, T, X);
  fclose(fp);

  /* Run the Baum-Welch algorithm and display the results */
  Baum_Welch(&dhmm, T, X, alpha, beta, gamma, xi);

  fp = fopen("t-a.out", "a");
  fprintf(fp, "\n* Baum-Welch algorithm:  maximum likelihood estimate:\n");
  write_DHMM(fp, &dhmm);
  fclose(fp);

  /* Run the penalized Baum-Welch algorithm and display the results */
  Penalized_Baum_Welch(&dhmm, &prior, T, X, alpha, beta, gamma, xi);

  fp = fopen("t-a.out", "a");
  fprintf(fp, "\n* Penalized Baum-Welch algorithm: posterior mode:\n");
  write_DHMM(fp, &dhmm);
  fclose(fp);

  /* Run the approximate Bayes algorithm and display the results */
  Approximate_Bayes(&dhmm, &prior, T, X, alpha, beta, gamma, xi);

  fp = fopen("t-a.out", "a");
  fprintf(fp, "\n* Approximate Bayes algorithm: approximate posterior distribution:\n");
  write_DHMM(fp, &dhmm);
  fclose(fp);


  /* Free the memory we used and close the filehandle */
  free_array_double_2d(alpha, 1, dhmm.I, 1, T);
  free_array_double_2d( beta, 1, dhmm.I, 1, T);
  free_array_double_2d(gamma, 1, dhmm.I, 1, T);
  free_array_double_3d(   xi, 1, dhmm.I, 1, dhmm.I, 1, T);
  free_array_int_1d(    X, 1, T);
  free_DHMM(&dhmm);
  free_DHMM(&model);
  free_DHMM(&prior);

}
