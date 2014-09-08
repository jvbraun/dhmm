/*      Program:  test-utilities.c
 *      Version:  0.1
 *         Date:  3 March 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Test the utility functions.
 */

static char rcsid[] = "$Id: tests.nw,v 1.2 1999/04/10 17:54:03 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdio.h>
#include "utilities.h"

/*********************************************************************
 * PROGRAM 
 ********************************************************************/

void main () {
  int i, j;  /* index variables */

  int *X;
  double *pi, **A, ***xi;

  double *p, *u;
  unsigned long seed;

  FILE *fp2;

  /* Set the random number generator seed        */ 
  set_seed(&seed, 0);

  fp2 = fopen("t-u.out", "w");
  fprintf(stdout, "Writing test results to 't-u.out'...\n");

  /* Test the uniform and gamma deviate generators */
  fprintf(fp2, "Testing uniform deviate generation, here are three...\n");
  for (i=1; i<=3; i++) {
    fprintf(fp2, "%lf ", uniform_deviate(&seed));
  }
  fprintf(fp2, "\n");

  fprintf(fp2, "Testing gamma deviate generation, here are three...\n");
  for (i=1; i<=3; i++) {
    fprintf(fp2, "%lf ", gamma_deviate(&seed, 100));
  }
  fprintf(fp2, "\n");

  /* Test the special functions */
  fprintf(fp2, "Testing digamma function...\n");
  fprintf(fp2, "digamma(1.0) = %lf  digamma(2.0) = %lf\n", digamma(1.0), digamma(2.0));

  /* Test array allocation, assignment, and deallocation */
  fprintf(fp2, "Testing array allocation...\n");

  X = array_int_1d(1,100);
  pi = array_double_1d(1,3);
  A = array_double_2d(1, 3, 1, 4);
  xi = array_double_3d(1, 3, 1, 4, 1, 100);

  fprintf(fp2, "Making some array assignments...\n");
  X[1] = 1;
  X[2] = 2;
  A[1][4] = 1.0;
  xi[2][3][99] = 1e-6;

  fprintf(fp2, "Testing array deallocation...\n");

  free_array_int_1d(X, 1, 100);
  free_array_double_1d(pi, 1, 3);
  free_array_double_2d(A, 1, 3, 1, 4);
  free_array_double_3d(xi,1, 3, 1, 4, 1, 100);

  /* Test the Dirichlet deviate generator */
  u = array_double_1d(1,5);
  p = array_double_1d(1,5);

  u[1] = 1; u[2] = 5; u[3] = 1; u[4] = 1; u[5] = 1;

  fprintf(fp2, "Testing Dirichlet deviate generation, here are three...\n");
  fprintf(fp2, "Prior u: ");
  for (j=1; j<=5; j++) {
    fprintf(fp2, "%lf ", u[j]);
  }
  fprintf(fp2, "\n");

  for (i=1; i<=3; i++) {
    fprintf(fp2, "Deviate: ");
    Dirichlet_deviate(&seed, u, 1, 5, p);
   for (j=1; j<=5; j++) {
      fprintf(fp2, "%lf ", p[j]);
        }
    fprintf(fp2, "\n");
  }
  fprintf(fp2, "\n");

  free_array_double_1d(u, 1, 5);
  free_array_double_1d(p, 1, 5);

  fprintf(fp2, "Looks like everything went okay...\n");

  fprintf(fp2, "End of tests.\n");
  fclose(fp2);

}
