/*      Program:  dhmm.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Provide utility functions for working with the 
 *                hidden Markov model structure.
 */

static char rcsid[] = "$Id: dhmm.nw,v 1.3 1999/05/04 21:56:37 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include "dhmm.h"
#include "utilities.h"

/*********************************************************************
 * FUNCTIONS 
 ********************************************************************/

/*           free_DHMM(&dhmm);
 * 
 * Purpose:  Free memory allocated to a DHMM structure.
 *
 *   To do:  1.  Check on deallocating the constants?
 *           2.  Provide some return for success.
 *
 *   Input:  pointer to DHMM structure that has been allocated
 *  Output:  none
 *  Return:  none
 */

void free_DHMM(DHMM *pdhmm) 
{
  free_array_double_2d(pdhmm->A, 1, pdhmm->I, 1, pdhmm->I);
  free_array_double_2d(pdhmm->B, 1, pdhmm->I, 1, pdhmm->M);
  free_array_double_1d(pdhmm->pi, 1, pdhmm->I);
}

/*           write_DHMM(stdout, &dhmm);
 * 
 * Purpose:  Write the DHMM structure to a filehandle in a uniform
 *           structure.
 *
 *   To do:  1.  Add return value.
 *
 *   Input:  pointer to output file, pointer to DHMM structure
 *  Output:  writes the information in the structure to the file 
 *  Return:  none
 */

void write_DHMM(FILE *fp, DHMM *pdhmm)
{
  int i, j, m;

  fprintf(fp, "I: %d\n", pdhmm->I);
  fprintf(fp, "M: %d\n", pdhmm->M);

  fprintf(fp, "A:\n");
  for (i=1; i<=pdhmm->I; i++) {
    for (j=1; j<=pdhmm->I; j++) {
      fprintf(fp, "%lf ", pdhmm->A[i][j]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "B:\n");
  for (i=1; i<=pdhmm->I; i++) {
    for (m=1; m<=pdhmm->M; m++) {
      fprintf(fp, "%lf ", pdhmm->B[i][m]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "pi:\n");
  for (i=1; i<=pdhmm->I; i++) {
    fprintf(fp, "%lf ", pdhmm->pi[i]);
  }
  fprintf(fp, "\n");  /* this trailing newline makes reading easier */
}

/*           read_DHMM(stdin, &dhmm);
 *
 * Purpose:  Read a DHMM structure from a file in the format 
 *           provided by write_DHMM.
 *
 *    Note:  This function allocates space for a DHMM structure.
 *
 *   Input:  pointer to input file, pointer to DHMM structure
 *  Output:  none
 *  Return:  none
 */

void read_DHMM(FILE *fp, DHMM *pdhmm) 
{
  int i, j, m, I, M;

  /* get the dimensions of the DHMM structure */
  fscanf(fp, "I: %d\n", &I);
  fscanf(fp, "M: %d\n", &M);

  /* allocate space */
  new_DHMM(pdhmm, I, M);

  /* continue to read in the information */
  fscanf(fp, "A:\n");
  for (i=1; i<=I; i++) {
    for (j=1; j<=I; j++) {
      fscanf(fp, "%lf", &(pdhmm->A[i][j]));  
    }
  fscanf(fp, "\n");
  }

  fscanf(fp, "B:\n");
  for (i=1; i<=I; i++) {
    for (m=1; m<=M; m++) {
      fscanf(fp, "%lf", &(pdhmm->B[i][m]));
    }
    fscanf(fp, "\n");
  }

  fscanf(fp, "pi:\n");
  for (i=1; i<=I; i++) {
    fscanf(fp, "%lf", &(pdhmm->pi[i]));
  }
}

/*           initialize_DHMM(&seed, &dhmm);
 *
 * Purpose:  Randomly initialize a DHMM structure
 *           with given number of states and observation types.
 *
 *    Note:  The DHMM must already be allocated.
 *  
 *   Input:  random number generator seed, 
 *           pointer to initialized DHMM structure, 
 *  Output:  none
 *  Return:  none
 */

void initialize_DHMM(unsigned long *seed, DHMM *pdhmm)
{
  int i, j, m, I, M;
  double sum;

  /* get the model dimensions */
  I = pdhmm->I;
  M = pdhmm->M;

  /* initialize pi, the initial state vector */
  sum = 0.0;
  for (i=1; i<=I; i++) {
    pdhmm->pi[i] = drand(seed);
    sum += pdhmm->pi[i];
    }
  for (i=1; i<=I; i++) {
    pdhmm->pi[i] /= sum;
  }

  /* initialize A, the transition matrix */
  for (i=1; i<=I; i++) {
    sum = 0.0;
    for (j=1; j<=I; j++) {
      pdhmm->A[i][j] = drand(seed);
      sum += pdhmm->A[i][j];
    }
    for (j=1; j<=I; j++) {
      pdhmm->A[i][j] /= sum;    /* normalize the rows */
    }
  }

  /* initialize B, the observation emission matrix */
  for (i=1; i<=I; i++) {
    sum = 0.0;
    for (m=1; m<=M; m++) {
      pdhmm->B[i][m] = drand(seed);
      sum += pdhmm->B[i][m];
    }
    for (m=1; m<=I; m++) {
      pdhmm->B[i][m] /= sum;    /* normalize the rows */
    }
  }
}

/*           DHMM_deviate(&seed, &prior, &dhmm);
 *
 * Purpose:  Generates realizations of discrete hidden Markov models
 *               given a prior distribution which consists of independent
 *               Dirichlet distributions on rows.
 *
 *    Note:  Both structures must be allocated.
 *
 *   Input:  random number generator seed, pointer to DHMM structure 
 *           with prior, pointer to DHMM structure for the realization
 *  Output:  none
 *  Return:  none
 */

void DHMM_deviate(unsigned long *seed, DHMM *pprior, DHMM *pdhmm)
{
  int i, I, M;

  /* get the model structure from the prior distribution */
  I = pprior->I;
  M = pprior->M;

  /* generate the rows of the realization 
   * each row of the realization is independent of the other rows,
   * and each row is generated as a realization from a Dirichlet
   * process with parameter vector given in the prior
   */
  
  /* pi */
  Dirichlet_deviate(seed, pprior->pi, 1, I, pdhmm->pi);

  /* A */
  for (i=1; i<=I; i++) {
    Dirichlet_deviate(seed, pprior->A[i], 1, I, pdhmm->A[i]);
  }

  /* B */
  for (i=1; i<=I; i++) {
    Dirichlet_deviate(seed, pprior->B[i], 1, M, pdhmm->B[i]);
  }
}

/*           copy_DHMM(&original, &copy);
 *
 * Purpose:  Copy the original DHMM structure into a copy DHMM structure.
 *
 *    Note:  It is the caller's responsibility to make sure
 *           that the size of the two DHMM's is the same.
 *
 *   Input:  pointer to two initialized DHMM structures
 *  Output:  none
 *  Return:  none
 */

void copy_DHMM(DHMM *poriginal, DHMM *pcopy)
{
  int i, j, m, I, M;

  /* get the model dimensions */
  I = poriginal->I;
  M = poriginal->M;

  /* copy pi */
  for (i=1; i<=I; i++) {
    pcopy->pi[i] = poriginal->pi[i];
  }

  /* copy A */
  for (i=1; i<=I; i++) {
    for (j=1; j<=I; j++) {
          pcopy->A[i][j] = poriginal->A[i][j];
        }
  }

  /* copy B */
  for (i=1; i<=I; i++) {
    for (m=1; m<=M; m++) {
          pcopy->B[i][m] = poriginal->B[i][m];
        }
  }
}

/*           new_DHMM(&dhmm, I, M);
 *
 * Purpose:  Allocate space for a DHMM with I states and with
 *           M observation types.
 *
 *   Input:  pointer to uninitialized DHMM structure
 *  Output:  none
 *  Return:  none
 */

void new_DHMM(DHMM *pdhmm, int I, int M)
{
  /* create the structure */
  pdhmm->I  = I;
  pdhmm->M  = M;

  pdhmm->A  = (double **) array_double_2d(1, I, 1, I);
  pdhmm->B  = (double **) array_double_2d(1, I, 1, M);
  pdhmm->pi = (double *) array_double_1d(1, I);
}


