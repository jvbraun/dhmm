/*      Program:  sequence.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Provide utility functions for working with the 
 *                discrete sequences.
 */

static char rcsid[] = "$Id: sequence.nw,v 1.3 1999/04/20 22:17:04 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "dhmm.h"
#include "utilities.h"


/*********************************************************************
 * FUNCTIONS 
 ********************************************************************/

/*           read_sequence(stdin, &T, &X);
 *
 * Purpose:  Read a sequence in the format of write_sequence from
 *           the file.  Set T to the sequence length, allocate space
 *           for the data, and read the data into X.
 *
 *   To do:  1.  Check whether X is already allocated?
 *
 *   Input:  pointer to input file, pointer to integer, 
 *           pointer to 1-d array pointer
 *  Output:  none
 *  Return:  none
 */

void read_sequence(FILE *fp, int *pT, int **pX)
{
  int t;
  int *X;

  fscanf(fp, "T: %d\n", pT);
  
  X = (int *) array_int_1d(1, *pT);

  fscanf(fp, "X:\n");
  for (t=1; t<=*pT; t++) {
    fscanf(fp, "%d", &X[t]);
  }
  fscanf(fp, "\n");     /* useful for batch training */ 
 
  *pX = X;              /* here is the data */
}

/*           write_sequence(stdout, T, X);
 *
 * Purpose:  Write sequence information to a filehandle in a uniform
 *           format.
 *
 *   Input:  pointer to output file, integer length of sequence, 
 *           pointer to sequence
 *  Output:  writes the information to the file
 *  Return:  none
 */

void write_sequence(FILE *fp, int T, int *X)
{
  int t;

  fprintf(fp, "T: %d\n", T);

  fprintf(fp, "X:\n");
  for (t=1; t<=T; t++) {
    fprintf(fp, "%d ", X[t]);
  }
  fprintf(fp, "\n");
}

/*           generate_sequence(&seed, dhmm, T, &X); 
 *
 * Purpose:  Allocate space for a sequence of length T, then fill T
 *           with a realization from the hidden Markov model contained
 *           in dhmm.
 *
 *   To do:  1.  Check whether X is already allocated?
 *
 *   Input:  pointer to DHMM structure, integer length of sequence, 
 *           pointer to data array
 *  Output:  none
 *  Return:  none
 */

void generate_sequence(unsigned long *seed, DHMM *pdhmm, int T, int **pX)
{
  int i, j, t, *X;

  X = (int *) array_int_1d(1, T);

  /* generate the initial state */
  i = multinomial_deviate(seed, pdhmm->pi, 1, pdhmm->I);

  /* generate the initial observation given we are in state i */
  X[1] = multinomial_deviate(seed, pdhmm->B[i], 1, pdhmm->M);

  /* now generate the rest of the sequence */
  for (t=2; t<=T; t++) {

    /* given we were in state i, which state shall we move to */
    j = multinomial_deviate(seed, pdhmm->A[i], 1, pdhmm->I);

    i = j; /* now state i represents the state we moved to */
  
    /* generate the observation given we are in state i */
    X[t] = multinomial_deviate(seed, pdhmm->B[i], 1, pdhmm->M);

  }
        
  /* here is the pointer to the data */
  *pX = X;
}


