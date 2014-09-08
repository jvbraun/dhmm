/*      Program:  utilities.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Provide general utility functions for use in 
 *                constructing discrete hidden Markov models.  
 *              
 *     Contents:  Utility error reporting subroutine
 *                Dynamic array memory freeing
 *                Dynamic array memory allocation
 *                Digamma function (psi function)
 *                Pseudorandom number generators
 */

static char rcsid[] = "$Id: utilities.nw,v 1.3 1999/07/24 07:05:47 Jerome Braun Exp $";

/*********************************************************************
 *  INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "utilities.h"

/*********************************************************************
 * FUNCTIONS 
 ********************************************************************/

/*           utilities_error(error_string);
 * 
 * Purpose:  Report error messages in array allocation.
 *
 *   Notes:  Inspired by
 *               "Numerical Recipes for C" by Press, Flannery, 
 *               Teukolsky, and Vetterling. 
 *
 *   Input:  an error message in a char array
 *  Output:  print the error message to stderr and exit with
 *           error level 1 
 *  Return:  none
 */

void utilities_error(char error_string[])
{
  void exit();
  fprintf(stderr, "error: %s\n", error_string);
  exit(1);
}

/*              parray = array_int_1d(i1, i2);
 *
 * Purpose:  Functions for dynamic allocation of memory 
 *           for multidimensional arrays.
 *
 *   Notes:  Inspired by comp.lang.c FAQ (Section 6.16) at
 *               http://www.eskimo.com/~scs/C-faq/top.html
 *           and 
 *               "Numerical Recipes for C" by Press, Flannery, 
 *               Teukolsky, and Vetterling. 
 * 
 *           I prefer the flexibility of using a non-zero lower index
 *           for the arrays.  It is the caller's responsibility 
 *           to obey the index restrictions.
 *
 *           Section 7.31 of the FAQ points out that calloc uses
 *           an all-bits-zero filling which does not necessarily
 *           guarantee floating point zero.
 *
 *   Input:  first and last array index for each dimension
 *  Output:  calls utilities_error if the malloc fails
 *  Return:  pointer to array with unitialized, allocated space
 */

int *array_int_1d (int i1, int i2)
{
  int *parray;

  parray = (int *) malloc( (i2-i1+1) * sizeof(int) );
  if (!parray) utilities_error("memory allocation failure in array_int_1d()");

  return parray - i1;
}

double *array_double_1d (int i1, int i2)
{
  double *parray;

  parray = (double *) malloc( (i2-i1+1) * sizeof(double) );
  if (!parray) utilities_error("memory allocation failure in array_double_1d()");

  return parray - i1;
}

double **array_double_2d (int i1, int i2, int j1, int j2)
{
  int i;
  double **parray;

  parray = (double **) malloc( (i2-i1+1) * sizeof(double *) );
  if (!parray) utilities_error("memory allocation failure (1) in array_double_2d()");
  parray -= i1;

  for (i=i1; i<=i2; i++) {
    parray[i] = (double *) malloc( (j2-j1+1) * sizeof(double));
    if (!parray[i]) utilities_error("memory allocation failure (2) in array_double_2d()");
    parray[i] -= j1;
  }
 
  return parray;
}

double ***array_double_3d (int i1, int i2, int j1, int j2, int k1, int k2)
{
  int i, j;
  double ***parray;

  parray = (double ***) malloc( (i2-i1+1) * sizeof(double **) );
  if (!parray) utilities_error("memory allocation failure (1) in array_double_3d()");
  parray -= i1;

  for (i=i1; i<=i2; i++) {
    parray[i] = (double **) malloc( (j2-j1+1) * sizeof(double *) );
    if (!parray[i]) utilities_error("memory allocation failure (2) in array_double_3d()");
    parray[i] -= j1;

    for (j=j1; j<=j2; j++) {
      parray[i][j] = (double *) malloc( (k2-k1+1) * sizeof(double) );
      if (!parray[i][j]) utilities_error("memory allocation failure (3) in array_double_3d()");
      parray[i][j] -= k1;
    }
  }
  
  return parray;
}

/*           free_array_int_1d(a, i1, i2);
 *
 * Purpose:  Frees memory allocated by the array allocation function
 *           array_int_1d.
 *
 *    Note:  It is the caller's responsibility to make sure the correct
 *           indices are specified.
 *
 *   Input:  array pointer and initial and last indices
 *  Output:  none
 *  Return:  none
 */

void free_array_int_1d(int *parray, int i1, int i2)
{
  free((char *) (parray + i1));
}

void free_array_double_1d(double *parray, int i1, int i2)
{
  free((char *) (parray + i1));
}

void free_array_double_2d(double **parray, int i1, int i2, int j1, int j2)
{
  int i;
  
  for (i=i2; i>=i1; i--) {
    free((char *) (parray[i] + j1));
  }
  free((char *) (parray + i1));
}

void free_array_double_3d(double ***parray, int i1, int i2, int j1, int j2, int k1, int k2)
{
  int i, j;

  for (i=i2; i>=i1; i--) {
    for (j=j2; j>=j1; j--) {
      free((char *) (parray[i][j] + k1));
    }
  }
  free((char *) (parray + i1));
}

/*           a = digamma(x);
 * 
 * Purpose:  Calculate the value the psi or digamma function, namely 
 *           the d/dx log Gamma (x).
 *
 *    Note:  From Algorithm AS 103 Applied Statistics (1976) Vol. 25, No. 3
 *
 *   To do:  1. Obtain a true double precision version of this function... 
 *           2. Add argument-checking.
 *
 *   Input:  positive double precision x
 *  Output:  none
 *  Return:  value of derivative of the log gamma function at x
 *
 */

double digamma (double x) 
{
  double result, y, r;

  /* Now a set of magic numbers from the Stirling expansion of the function */
  double s = 1.0e-05;
  double c = 8.5;
  double s3 = 8.333333333e-02;
  double s4 = 8.333333333e-03;
  double s5 = 3.968253968e-03;
  double d1 = -0.5772156649;

  y = x;
  result = 0.0;
  /* Return an error if the argument is non-positive... */

  /* use approximation if argument is <= s */
  if (y <= s) {
    result = d1 - 1.0/y;
        return result;
  }

  /* reduce to digamma(x+n) where (x+n) >= c */
  /* do until y >= c */
  while (y < c) {
    result = result - 1.0/y;
        y += 1.0;
  }

  /* use de Moivre's expansion if argument >= c */
  r = 1.0/y;
  result = result + log(y) - 0.5 * r;
  r = r * r;
  result = result - r * (s3 - r * (s4 - r*s5));
  return result;
}

/*           u = uniform_deviate(&seed);
 *
 * Purpose:  Generate pseudorandom numbers from the uniform (0,1) 
 *           distribution.
 *
 *    Note:  This is just a wrapper function to ease changing the
 *           underlying generator. 
 *
 *   Input:  pointer to unsigned long seed
 *  Output:  none
 *  Return:  pseudorandom deviate between 0 and 1
 */

double uniform_deviate(unsigned long *seed) 
{
   return drand(seed);
}

/*           set_seed(&seed, k);
 * 
 * Purpose:  Set the random number seed by either reading a 
 *           .SEED file if it exists, or by taking the input
 *           integer and transforming it into a seed value.
 *
 *    Note:  This is a wrapper function to make changes easier.
 *           k must be a value such that 0 < 4*k+1 < 2**32.
 *
 *   To do:  1.  Add argument checking.
 *
 *   Input:  pointer to unsigned long seed, unsigned long
 *  Output:  none
 *  Return:  none
 */  

void set_seed(unsigned long *seed, unsigned long k) 
{
  double u;
  FILE *SEED;

  /* if a .SEED file does not exist, set the seed per the specifications
   * for drand() and initialize the generator 
   */
  if ( (SEED = fopen(".SEED", "r")) == NULL ) {
    *seed = 4 * k + 1;
    u = uniform_deviate(seed);
  }

  /* otherwise, read in the seed value */
  else {
    fscanf(SEED, "%lu", seed);
    fclose(SEED);
  }

}

/*           save_seed(&seed);
 * 
 * Purpose:  Save the random number seed to a .SEED file.
 *
 *   To do:  1.  Add argument checking.
 *
 *   Input:  pointer to unsigned long seed
 *  Output:  writes seed value to .SEED
 *  Return:  none
 */  

void save_seed(unsigned long *seed) 
{
  FILE *SEED;

  /* open the .SEED file and write out the seed */
  SEED = fopen(".SEED", "w");
  fprintf(SEED, "%lu", *seed);
  fclose(SEED);
}

/*           u = multinomial_deviate(&seed, &p, first, last);
 *
 * Purpose:  Generate pseudorandom multinomials from the probability
 *           distribution given by the input vector, with indices
 *           ranging from low to high.
 *
 *    Note:  The caller is responsible for making sure that the
 *           array contains a discrete probability distribution.
 *
 *   Input:  pointer to unsigned long seed, pointer to an array
 *           of probabilities, first index of the array, last index
 *           of the array
 *  Output:  none
 *  Return:  pseudorandom deviate between start and finish inclusive
 *           representing the multinomial class chosen
 */

int multinomial_deviate(unsigned long *seed, double *p, int first, int last) 
{
  int i;
  double r;

  /* default to the first index */
  i = first;
  r = uniform_deviate(seed);

  /* the following code is equivalent to laying the categories out as
   * bins on the unit interval, then determining which bin r fell into
   */
  while (r >= p[i]) {
    r -= p[i];
    i++;

    /* do a sanity check */
    if (i > last) utilities_error("exceeded array bound in mutinomial_deviate()");
  }

  return i;

}

/*           g = gamma_deviate(&seed, a);
 *
 * Purpose:  Generate pseudorandom numbers from the gamma(a) 
 *           distribution.
 *
 *    Note:  This is just a wrapper function to ease changing the
 *           underlying generator.
 *
 *   Input:  pointer to unsigned long seed, double shape parameter a
 *  Output:  none
 *  Return:  double gamma deviate
 */

double gamma_deviate(unsigned long *seed, double a)
{
  return gll(seed, a);
}

/*           Dirichlet_deviate(&seed, u, i1, i2, p);
 *
 * Purpose:  Generates psuedorandom vectors from the Dirichlet(u) 
 *           distribution.
 * 
 * Algorithm:  This algorithm uses the standard result that a vector
 *           of independent gamma deviates with shape parameters 
 *           given by u, normalized by their sum, is distributed
 *           as Dirichlet(u).
 *
 *   Input:  pointer to unsigned long seed, pointer to parameter array,
 *           integer first and last index, pointer to result array
 *  Output:  replaces the result array with the generated deviate
 *  Return:  none
 */

void Dirichlet_deviate(unsigned long *seed, double *u, int i1, int i2, double *p)
{
  int i;      /* index variable */
  double sum; /* the sum of the individual gammas */

  /* generate a set of independent gamma distributed random
   * variates with the same scale (default here) and with
   * shape parameter given by the values in u
   */

  sum = 0.0;
  for (i=i1; i<=i2; i++) {
    p[i] = gamma_deviate(seed, u[i]);
        sum += p[i];
  }

  /* then normalize the values of the vector by their sum */
  for (i=i1; i<=i2; i++) {
    p[i] /= sum;
  }
}

/*********************************************************************
 * LIBRARY FUNCTIONS
 * 
 * Source:  Win_rand 1.0
 * Author:  Ernst Stadlober
 *    URL:  http://random.mat.sbg.ac.at/ftp/pub/data/winrand.zip
 * 
 *   Note:  The following functions were drawn from the Win_rand 1.0
 *          library of random number generators.  The library contains 
 *          many implementations of Luc Devroye's pseudorandom number 
 *          generators.This library also contains a nice GUI front end 
 *          for exploring the properties of the generators.  
 ********************************************************************/

/******************************************************************
 *                                                                *
 *                    (0,1) Uniform Distribution                  *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION: - drand   samples a random number from the (0,1)     *
 *             uniform distribution by a multiplicative           *
 *             congruential method of the form                    *
 *             z(i+1) = a * z(i) (mod m)                          *
 *             a=663608941 (=0X278DDE6DL), m=2**32.               *
 *           - Before the first call of the generator an integer  *
 *             value z(0)=4*k+1 where 0 < 4*k+1 < m , has to be   *
 *             initialized. ( e.g  *seed = 1; u=drand(seed);)   *
 *             Then z(i) assumes all different values 0<4*k+1<m   *
 *             during a full period 2**30 of drand                *
 ******************************************************************/

double drand(unsigned long *seed)
{
 *seed *= 0X278DDE6DL;                  /* z(i+1)=a*z(i) (mod 2**32) */
 return(2.32830643653870e-10 * *seed);  /* u(i) = z(i)*2**(-32) */
 }


/******************************************************************
 *                                                                *
 *   Gamma Distribution - Acceptance Rejection with log-logistic  *
 *                        envelopes                               *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION :   - gll produces a sample from the standard gamma   *
 *                distribution with parameter  a > 0.             *
 * REFERENCE :  - R.C.H. Cheng (1977): The generation of gamma    *
 *                variables with non-integral shape parameter,    *
 *                Appl. Statist. 26, 71-75.                       *
 * SUBPROGRAM : - drand(seed) ... (0,1)-Uniform generator with    *
 *                unsigned long integer *seed.                    *
 *                                                                *
 * Implemented by R. Kremer, 1990                                 *
 ******************************************************************/

double gll(unsigned long *seed, double a)
{
 static double aa,bb,cc,a_in = -1.0;
 double u1,u2,v,r,z,gl;

 if (a != a_in)
   {
    a_in = a;
    aa = (a > 1.0)? sqrt(a + a - 1.0) : a;
    bb = a - 1.386294361;
    cc = a + aa;
   }
 for(;;)
   {
    u1 = drand(seed);
    u2 = drand(seed);
    v = log(u1 / (1.0 - u1)) / aa;
    gl = a * exp(v);
    r = bb + cc * v - gl;
    z = u1 * u1 * u2;
    if (r + 2.504077397 >= 4.5 * z) break;
    if (r >= log(z)) break;
   }
 return(gl);
}


