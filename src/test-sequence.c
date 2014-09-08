/*      Program:  test-sequence.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Test the DHMM and sequence utility functions.
 */

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include "dhmm.h"
#include "utilities.h"

/*********************************************************************
 * PROGRAM 
 ********************************************************************/

void main (void) {
  DHMM dhmm, dhmm2, dhmm3;
  int I, M, T, status, i, j, m, t;
  int *X, *X2;

  FILE *fp, *fp2;
  unsigned long seed;  /* pseudorandom generator seed */

  /* set the random number generator seed        */ 
  set_seed(&seed, 0);

  I = 2;               /* number of states */
  M = 2;               /* number of observation types */
  T = 20;              /* sequence length */

  fprintf(stdout, "Writing test results to 't-s.out'...\n");

  fp2 = fopen("t-s.out", "w");

  /* allocate space for the DHMM's */
  fprintf(fp2, "Testing allocating space for DHMM structures...\n");
  new_DHMM(&dhmm, I, M);
  new_DHMM(&dhmm3, I, M);

  /* Test the DHMM structure utilities */
  fprintf(fp2, "Testing initializing a discrete hidden Markov model...\n");
  initialize_DHMM(&seed, &dhmm);

  fprintf(fp2, "Trying the DHMM copying utility...\n");
  copy_DHMM(&dhmm, &dhmm3);

  fprintf(fp2, "Testing working with the DHMM structure...\n");
  dhmm.A[1][1] = 0.5; dhmm.A[1][2] = 0.5;
  dhmm.A[2][1] = 0.1; dhmm.A[2][2] = 0.9;

  dhmm.B[1][1] = 0.7; dhmm.B[1][2] = 0.3;
  dhmm.B[2][1] = 0.2; dhmm.B[2][2] = 0.8;

  dhmm.pi[1] = 0.4; dhmm.pi[2] = 0.6;

  fprintf(fp2, "Testing writing the DHMM structure to 'd-s.dhm'...\n");
  fp = fopen("t-s.dhm", "w");
  write_DHMM(fp, &dhmm);
  fclose(fp);

  fprintf(fp2, "Testing reading the DHMM structure from 'd-s.dhm'...\n");
  fp = fopen("t-s.dhm", "r");
  read_DHMM(fp, &dhmm2);   /* allocates space for dhmm2 */
  fclose(fp);

  status = 0;
  for (i=1; i<=I; i++) {
  if (dhmm.pi[i] != dhmm2.pi[i]) status = 1;
    for (j=1; j<=I; j++) {
      if (dhmm.A[i][j] != dhmm2.A[i][j]) status = 1;
    }
    for (m=1; m<=M; m++) {
      if (dhmm.B[i][m] != dhmm2.B[i][m]) status = 1;
    }
  }

  fprintf(fp2, "Read and write of DHMM structure ");
  if (status == 0) {
    fprintf(fp2, "succeeded...\n");
  }
  else {
    fprintf(fp2, "failed...\n");
  }
  
  /* Test the sequence utilities */
  fprintf(fp2, "Testing generating a sequence...\n");
  generate_sequence(&seed, &dhmm, T, &X);

  fprintf(fp2, "Testing writing out the sequence to 't-s.seq'...\n");
  fp = fopen("t-s.seq", "w");
  write_sequence(fp, T, X);
  fclose(fp);

  fprintf(fp2, "Testing reading the sequence from 't-s.seq'...\n");
  fp = fopen("t-s.seq", "r");
  read_sequence(fp, &T, &X2);
  fclose(fp);

  status = 0;
  for (t=1; t<=T; t++) if (X[t] != X2[t]) status = 1;

  fprintf(fp2, "Read and write of sequence ");
  if (status == 0) {
    fprintf(fp2, "succeeded...\n");
  }
  else {
    fprintf(fp2, "failed...\n");
  }

  /* Test the pseudorandom DHMM generator */
  fprintf(fp2, "Testing generating a DHMM...\n");
  DHMM_deviate(&seed, &dhmm, &dhmm3);

  fprintf(fp2, "End of tests.\n");
  fclose(fp2);
}
