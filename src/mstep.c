/*      Program:  mstep.c
 *      Version:  0.1
 *         Date:  1 February 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Routines to carry out the maximization steps
 *                for discrete hidden Markov models.
 */

static char rcsid[] = "$Id: mstep.nw,v 1.3 1999/07/10 09:05:17 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES
 ********************************************************************/

#include <stdio.h>      /* else the function prototypes in dhmm.h fail */
#include <math.h>       /* need exp */
#include "dhmm.h"
#include "utilities.h"  /* need the digamma function */

/*********************************************************************
 * FUNCTIONS
 ********************************************************************/

/*           Classical_M_step(dhmm, epsilon, T, X, gamma, xi); 
 *
 * Purpose:  Perform the classical maximization step in the 
 *           Baum-Welch algorithm.
 *
 *    Note:  Changes dhmm in place.
 *
 *   To do:  1.  Amplify the above.
 *
 *   Input:  pointer to DHMM structure, double flattening constant,
 *           integer sequence length, data vector,
 *           gamma and xi arrays
 *  Output:  none
 *  Return:  none
 */

void Classical_M_step (DHMM *pdhmm, double epsilon, int T, int *X, double **gamma, double ***xi)
{
  int i, j, m, t;
  double sum, gammasum, numerator;

  /* reestimate pi */
  sum = 0.0;                                                     
  for (i=1; i <= pdhmm->I; i++) {
    pdhmm->pi[i] = epsilon + gamma[i][1];  /* initialize bin */
    sum += pdhmm->pi[i];
  }
  for (i=1; i <= pdhmm->I; i++) {
    pdhmm->pi[i] /= sum;                /* normalize */
  }  

  /* reestimate A */
  for (i=1; i <= pdhmm->I; i++) {

    /* get the denominator for the re-estimate of A */
    gammasum = 0.0;                    
    for (t=1; t<=T-1; t++)
      gammasum += gamma[i][t];

    /* now get the numerator and normalizing factor */
    sum = 0.0;
    for (j=1; j <= pdhmm->I; j++) {
      numerator = 0.0;             /* initialize the bin */
      for (t=1; t<=T-1; t++) {
        numerator += xi[i][j][t];
      }
      pdhmm->A[i][j] = epsilon + numerator/gammasum;
      sum += pdhmm->A[i][j];
    }

    for (j=1; j <= pdhmm->I; j++) {
      pdhmm->A[i][j] /= sum;            /* normalize */
    }

  }

  /* reestimate B */
  /* it is a little tricky in that we must use for the numerator */
  /* only those m which actually appeared                        */

  for (i=1; i <= pdhmm->I; i++) {
    for (m=1; m <= pdhmm->M; m++) {
      pdhmm->B[i][m] = epsilon;         /* initialize the bin */
    }
  }
  for (i=1; i <= pdhmm->I; i++) {       /* only those who got an observation get something */
    for (t=1; t<=T; t++) {
      pdhmm->B[i][X[t]] += gamma[i][t];
    }
  }
  for (i=1; i <= pdhmm->I; i++) {               /* normalize */
    sum = 0.0;
    for (t=1; t<=T; t++) {
      sum += gamma[i][t];
    }
    for (m=1; m <= pdhmm->M; m++) {
      pdhmm->B[i][m] /= sum; 
    }   
  }     
}

/*           Penalized_M_step(dhmm, prior, T, X, gamma, xi); 
 *
 * Purpose:  Perform the classical maximization step, but initializing
 *           the "bins" with the values in the prior.  This is equivalent
 *           to maximum a posterior inference with the softmax basis
 *           (MacKay).
 *
 *    Note:  Changes dhmm in place.
 *
 *   To do:  1.  Double-check this algorithm.
 *
 *   Input:  pointers to DHMM structures for model and priors, 
 *           double flattening constante,
 *           integer sequence length, data vector,
 *           gamma and xi arrays
 *  Output:  none
 *  Return:  none
 */

void Penalized_M_step (DHMM *pdhmm, DHMM *pprior, int T, int *X, double **gamma, double ***xi)
{
  int i, j, m, t;
  double sum, gammasum, numerator;

  /* reestimate pi */
  gammasum = 0.0;                                                     
  for (i=1; i <= pdhmm->I; i++) {
    gammasum += gamma[i][1];
  }

  /* move the estimate toward the prior */
  for (i=1; i <= pdhmm->I; i++) {
    pdhmm->pi[i] = pprior->pi[i] + gamma[i][1]/gammasum;
  }

  /* normalize */
  sum = 0.0;
  for (i=1; i <= pdhmm->I; i++) {
    sum += pdhmm->pi[i];
  }  

  for (i=1; i <= pdhmm->I; i++) {
    pdhmm->pi[i] /= sum;
  }

  /* reestimate A */
  for (i=1; i <= pdhmm->I; i++) {

    /* get the denominator for the re-estimate of A */
    gammasum = 0.0;
    for (t=1; t<=T-1; t++) {
      gammasum += gamma[i][t];
    }

    /* now get the numerator and the normalizing factor */
    sum = 0.0;                    
    for (j=1; j<=pdhmm->I; j++) {
      numerator = 0.0;                  /* initialize the bin */
      for (t=1; t<=T-1; t++) {
        numerator += xi[i][j][t];  /* sum over t */
      }
      pdhmm->A[i][j] = pprior->A[i][j] + numerator/gammasum;
      sum += pdhmm->A[i][j];
    }

    /* normalize the rows of A */
    for (j=1; j <= pdhmm->I; j++) {
      pdhmm->A[i][j] /= sum;
    }

  }

  /* reestimate B */
  /* it is a little tricky in that we must use for the numerator */
  /* only those m which actually appeared                        */

  for (i=1; i <= pdhmm->I; i++) {
    for (m=1; m <= pdhmm->M; m++) {
      pdhmm->B[i][m] = 0;
    }
  }

  /* only those who got an observation get something */
  for (i=1; i <= pdhmm->I; i++) {

    gammasum = 0.0;
    for (t=1; t<=T; t++) {
      pdhmm->B[i][X[t]] += gamma[i][t];
      gammasum += gamma[i][t];
    }

    /* now move the estimate toward the prior */
    for (m=1; m<=pdhmm->M; m++) {
      pdhmm->B[i][m] = pprior->B[i][m] + pdhmm->B[i][m]/gammasum;
    }

  }

  /* normalize the rows of B */
  for (i=1; i <= pdhmm->I; i++) {

    sum = 0.0;
    for (m=1; m <= pdhmm->M; m++) {
      sum += pdhmm->B[i][m];
    }

    /* calculate the normalized estimate */
    for (m=1; m <= pdhmm->M; m++) {
      pdhmm->B[i][m] /= sum; 
    }   

  }     
}

/*           Approximate_Bayes_M_step(posterior, prior, T, X, gamma, xi); 
 *
 * Purpose:  Perform the maximization step to minimize the free energy
 *           given alpha and beta (MacKay).
 *
 *    Note:  Changes posterior in place.
 *
 *   To do:  1.  Amplify the above.
 *
 *   Input:  pointers to DHMM structures for model and prior, 
 *           double flattening constante,
 *           integer sequence length, data vector,
 *           gamma and xi arrays
 *  Output:  none
 *  Return:  none
 */

void Approximate_Bayes_M_step (DHMM *pposterior, DHMM *pprior, int T, int *X, double **gamma, double ***xi)
{
  int i, j, m, t, I, M;

  /* get the model dimensions */
  I = pprior->I;
  M = pprior->M;

  /* calculate pi^* */
  for (i=1; i<=I; i++) {
    pposterior->pi[i] = pprior->pi[i] + gamma[i][1];    /* initialize bin with the prior value */
  }


  /* calculate A^* */
  for (i=1; i<=I; i++) {
    for (j=1; j<=I; j++) {
      pposterior->A[i][j] = pprior->A[i][j];    /* initialize the bin with the prior value */
      for (t=1; t<=T-1; t++) {
        pposterior->A[i][j] += xi[i][j][t];             /* sum over t */
      }
    }
  }

  /* B^* */
  /* it is a little tricky in that we must use for the numerator */
  /* only those m which actually appeared                        */

  for (i=1; i<=I; i++) {                
    for (m=1; m<=M; m++) {
      pposterior->B[i][m] = pprior->B[i][m];    /* initialize the bin with the prior value */
    }
  }
  for (i=1; i<=I; i++) {        /* only those who got an observation get something */
    for (t=1; t<=T; t++) {
      pposterior->B[i][X[t]] += gamma[i][t];
    }
  }

}

/*           Calculate_star(&posterior, &star);
 *
 * Purpose:  Perform the final part of the calculation for 
 *           the algorithm.
 *
 *   Input:  pointer to DHMM structure
 *  Output:  none
 *  Return:  none
 */

void Calculate_star(DHMM *pposterior, DHMM *pstar)
{
  int i, j, m, I, M;
  double sum;

  /* get the DHMM dimensions */
  I = pposterior->I;
  M = pposterior->M;

  /* subnormalize pi */
  sum = 0.0;
  for (i=1; i<=I; i++) {
    sum += pposterior->pi[i];
  }
  for (i=1; i<=I; i++) {
    pstar->pi[i] = exp( digamma(pposterior->pi[i]) - digamma(sum) );
  }

  /* subnormalize A */
  for (i=1; i<=I; i++) {
    sum = 0.0;                          /* accumulate the Dirichlet parameter */  
    for (j=1; j<=I; j++) {
      sum += pposterior->A[i][j];
    }
    for (j=1; j<=I; j++) {              /* now complete the calculation */
      pstar->A[i][j] = exp( digamma(pposterior->A[i][j]) - digamma(sum) ); 
    }
  }

  /* subnormalize B */
  for (i=1; i<=I; i++) {
    sum = 0.0;                  /* accumulate the Dirichlet parameter */
    for (m=1; m<=M; m++) {
      sum += pposterior->B[i][m];
    }
    for (m=1; m<=M; m++) {      /* now complete the calculation */
      pstar->B[i][m] = exp( digamma(pposterior->B[i][m]) - digamma(sum) );
    }
  }
}

