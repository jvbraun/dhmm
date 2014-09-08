/*      Program:  estep.c
 *      Version:  0.1
 *         Date:  1 February 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Routines to carry out the expectation steps
 *                for discrete hidden Markov models.
 */

static char rcsid[] = "$Id: estep.nw,v 1.2 1999/04/10 17:54:03 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES
 ********************************************************************/

#include <stdio.h>      /* else the function prototypes in dhmm.h fail */
#include <math.h>       /* needed for the log function */
#include "dhmm.h" 

/*********************************************************************
 * FUNCTIONS
 ********************************************************************/

/* The following functions are used to carry out the classical
 * Baum-Welch algorithm.  For practical use, it is necessary to
 * use the calculation with scaling, since even moderate sequence
 * lengths result in extremely tiny conditional probabilities.
 */

/*           Forward(dhmm, T, X, alpha);
 *
 * Purpose:  Calculate the forward probabilities in the Baum-Welch
 *           algorithm.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           array to hold the probabilities
 *  Output:  none
 *  Return:  none 
 */

void Forward(DHMM *pdhmm, int T, int *X, double **alpha)
{
  int i, j, t;
  double sum;

  /* initialize the recursion */
  for (i=1; i<=pdhmm->I; i++) {
    alpha[i][1] = pdhmm->pi[i] * pdhmm->B[i][X[1]];
  }

  /* carry out the recursion */
  for (t=2; t<=T; t++) {
    for (j=1; j<=pdhmm->I; j++) {  
      sum = 0.0;      
      for (i=1; i<=pdhmm->I; i++) {
        sum += alpha[i][t-1] * pdhmm->A[i][j];
      }
      alpha[j][t] = sum * pdhmm->B[j][X[t]];     /* safer to accumulate then multiply */
    }
  }

}

/*           Forward_Scaled(dhmm, T, X, alpha, scale);
 *
 * Purpose:  Calculate the forward probabilities in the Baum-Welch
 *           algorithm, using scaling to prevent underflow.
 *
 *  Method:  The scaling used here is the usual one which normalizes
 *           the alpha_{ij}^t for given t.  
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           array to hold the probabilities, array to hold scale
 *  Output:  none
 *  Return:  none 
 */

void Forward_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, double *scale)
{
  int i, j, t, I, M;
  double sum;

  /* get the model dimensions */
  I = pdhmm->I;
  M = pdhmm->M;

  /* initialize the recursion */
  scale[1] = 0.0;
  for (i=1; i<=I; i++) {
    alpha[i][1] = pdhmm->pi[i] * pdhmm->B[i][X[1]];
    scale[1] += alpha[i][1];
  }
  for (i=1; i<=I; i++) {
    alpha[i][1] /= scale[1];    /* rescale the alpha's */
  }

  for (t=2; t<=T; t++) {
    scale[t] = 0.0;
    for (j=1; j<=I; j++) {  
      sum = 0.0;      
      for (i=1; i<=I; i++) {
            sum += alpha[i][t-1] * pdhmm->A[i][j];
      }
      alpha[j][t] = sum * pdhmm->B[j][X[t]];     /* better to accumulate then multiply */
      scale[t] += alpha[j][t];
    }
    /* rescale the alpha's for each t */
        for (j=1; j<=I; j++) {
      alpha[j][t] /= scale[t];  
    }
  }

}

/*           Backward(dhmm, T, X, beta);
 *
 * Purpose:  Calculate the backward probabilities in the Baum-Welch
 *           algorithm.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           array to hold the probabilities
 *  Output:  none
 *  Return:  none 
 */

void Backward(DHMM *pdhmm, int T, int *X, double **beta)
{
  int i, j, t; 
  double sum;

  /* begin the recursion */
  for (i=1; i<=pdhmm->I; i++) {
    beta[i][T] = 1.0;
  }

  for (t=T-1; t>=1; t--) {
    for (i=1; i<=pdhmm->I; i++) {
      sum = 0.0;
      for (j=1; j<=pdhmm->I; j++) {
        sum += pdhmm->A[i][j] * beta[j][t+1] * pdhmm->B[j][X[t+1]];
      }
    beta[i][t] = sum;
    }
  }
}

/*           Backward_Scaled(dhmm, T, X, beta, scale);
 *
 * Purpose:  Calculate the backward probabilities in the Baum-Welch
 *           algorithm.
 *
 *    Note:  The scaling is assumed already computed and is passed
 *           into the function.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           array to hold the probabilities
 *  Output:  none
 *  Return:  none 
 */

void Backward_Scaled(DHMM *pdhmm, int T, int *X, double **beta, double *scale)
{
  int i, j, t, I, M; 
  double sum;

  /* get the model dimensions */
  I = pdhmm->I;
  M = pdhmm->M;

  /* begin the recursion */
  for (i=1; i<=I; i++) {
    beta[i][T] = 1.0 / scale[T];
  }

  for (t=T-1; t>=1; t--) {
    for (i=1; i<=I; i++) {
      sum = 0.0;
      for (j=1; j<=I; j++) {
        sum += pdhmm->A[i][j] * beta[j][t+1] * pdhmm->B[j][X[t+1]];
      }
    beta[i][t] = sum / scale[t];
    }
  }
}

/*           Normalize_xi(dhmm, T, X, alpha, beta, xi);
 *
 * Purpose:  Normalize the elements of the xi array.
 * 
 *   To do:  1.  Amplify on the above.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           forward and backward probability arrays,
 *           array for results
 *  Output:  none
 *  Return:  none 
 */

void Normalize_xi(double ***xi, int I, int M, int T)
{
  int i, j, t;
  double sum;

  for (t=1; t<=T-1; t++) {
    /* sum over i, j for each t */
        sum = 0.0;
    for (i=1; i<=I; i++) {
      for (j=1; j<=I; j++) {
        sum += xi[i][j][t];
          }
    }
    /* normalize them for each t */
    for (i=1; i<=I; i++) {
      for (j=1; j<=I; j++) {
        xi[i][j][t] /= sum;
      }
    }
  }
}

/*           Normalize_gamma(dhmm, T, X, alpha, beta, gamma);
 *
 * Purpose:  Normalize the gamma array.
 * 
 *   To do:  1.  Amplify on the above.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           forward and backward probability arrays,
 *           array for results
 *  Output:  none
 *  Return:  none 
 */

void Normalize_gamma(double **gamma, int I, int T)
{
  int i, t;
  double sum;

  for (t=1; t<=T; t++) {
    /* sum over the i for each t */
        sum = 0.0;
    for (i=1; i<=I; i++) {
      sum += gamma[i][t];
    }
    /* normalize them for each t */
    for (i=1; i<=I; i++) {
      gamma[i][t] /= sum;
    }
  }
}

/*           Accumulate_xi(&dhmm, T, X, alpha, beta, xi);
 *
 * Purpose:  Accumulate the elements of the xi array without normalizing
 *           them.
 * 
 *   To do:  1.  Amplify on the above.
 * 
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           forward and backward probability arrays,
 *           array for results, array for scaling
 *  Output:  none
 *  Return:  none 
 */

void Accumulate_xi(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double ***xi)
{
  int i, j, t, I;

  /* get model dimensions */
  I = pdhmm->I;

  /* calculate the probabilities of transitions */
  /* probability of transition from i->j at step t given X */
  for (t=1; t<=T-1; t++) {
    for (i=1; i<=I; i++) {
      for (j=1; j<=I; j++) {
        xi[i][j][t] = alpha[i][t] * pdhmm->A[i][j] * 
                              beta[j][t+1] * pdhmm->B[j][X[t+1]];
      }
    }
  }

}

/*           Accumulate_xi_Scaled(&dhmm, T, X, alpha, beta, xi, scale);
 *
 * Purpose:  Accumulate the elements of the xi array without normalizing
 *           them -- scaled version.
 * 
 *   To do:  1.  Amplify on the above.
 * 
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           forward and backward probability arrays,
 *           array for results, array for scaling
 *  Output:  none
 *  Return:  none 
 */

void Accumulate_xi_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double ***xi, double *scale)
{
  int i, j, t, I;

  /* get model dimensions */
  I = pdhmm->I;

  /* calculate the probabilities of transitions */
  /* probability of transition from i->j at step t given X */
  for (t=1; t<=T-1; t++) {
    for (i=1; i<=I; i++) {
      for (j=1; j<=I; j++) {
        xi[i][j][t] = alpha[i][t] * pdhmm->A[i][j] * 
                              beta[j][t+1] * pdhmm->B[j][X[t+1]] / scale[t+1];
          /* note the scaling factor */
      }
    }
  }
}

/*           Accumulate_gamma(&pdhmm, T, X, alpha, beta, gamma);
 *
 * Purpose:  Accumulate the values of the gamma array without normalizing
 *           them.
 *
 *   To do:  1.  Amplify on the above.
 *
 *   Input:  pointer to DHMM structure, 
 *           integer sequence length, sequence data,
 *           forward and backward probability arrays,
 *           array for results, array for scale
 *  Output:  none
 *  Return:  none 
 */

void Accumulate_gamma(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double **gamma)
{
  int i, t, I;

  /* get model dimension */
  I = pdhmm->I;

  /* calculate the probabilities of being in state i at t */
  for (t=1; t<=T; t++) {
    for (i=1; i<=I; i++) {
      gamma[i][t] = alpha[i][t] * beta[i][t];
    }
  }
}

/*           Log-Likelihood(alpha, I, T);
 *
 * Purpose:  Calculate the log-likelihood of the data assuming
 *           that the unscaled forward probabilities reside in 
 *           alpha.
 * 
 *   To do:  1.  Amplify on the above.
 *
 *   Input:  pointer to alpha array, integer number of states,
 *           integer number of observations
 *  Output:  none
 *  Return:  double precision log-likelihood
 */

double Log_Likelihood(double **alpha, int I, int T)
{
  int i;
  double LL;

  /* log-likelihood of the data given the current parameter estimates */
  LL = 0.0;
  for (i=1; i<=I; i++) {
    LL += alpha[i][T];
  }
  return log(LL);
}

/*           Log-Likelihood_Scaled(scale, T);
 *
 * Purpose:  Calculate the log-likelihood of the data assuming
 *           that the scaling of Forward_Scaled was used.
 *
 *  Method:  This scaling has the nice side effect that the
 *           sum is obtained for each t over given i and j.
 * 
 *   To do:  1.  Amplify on the above.
 *
 *   Input:  pointer to scale array, integer number of observations
 *  Output:  none
 *  Return:  double precision log-likelihood
 */

double Log_Likelihood_Scaled(double *scale, int T)
{
  int t;
  double LL;

  /* this scaling is also convenient, since it gives us automatically */
  /* the sum of the direct probabilities */
  LL = 0.0;
  for (t=1; t<=T; t++) {
    LL += log(scale[t]);
  }

  return LL;
}

/*           FE = Free_Energy(dhmm, prior, X, T, gamma, xi);
 *
 * Purpose:  Calculate the log free energy, int_{Q} log P(X S theta | U),
 *           between the approximation to the posterior and the true posterior.
 *
 *  Method:  The varying part of the free energy is calculated.
 *           We will use the fact that we have already calculated exp(E[a_ij]),
 *           exp(E[b_ij]), and exp(E[pi_i]) as part of the star calculation.
 *
 *   Input:  pointers to DHMM structures for model and prior,
 *           pointer to the free energy,
 *           data array and data length,
 *           gamma and xi arrays
 *  Output:  none
 *  Return:  double precision free energy
 */

double Free_Energy(DHMM *pstar, DHMM *pprior, int *X, int T, double **gamma, double ***xi)
{
  int i, j, m, t;
  int I, M;

  double FE;

  I = pstar->I;    /* number of states */
  M = pstar->M;    /* number of observation types */
 
  FE = 0.0;

  for (i=1; i<=I; i++) {
          FE += ( pprior->pi[i] - 1 ) * log(pstar->pi[i]);
  }

  for (i=1; i<=I; i++) {
    for (j=1; j<=I; j++) {
      FE += ( pprior->A[i][j] - 1 ) * log(pstar->A[i][j]);
    }
  }

  for (i=1; i <= pstar->I; i++) {
    for (m=1; m <= pstar->M; m++) {
      FE += ( pprior->B[i][m] - 1 ) * log(pstar->B[i][m]);
    }
  }

  for (i=1; i <= pstar->I; i++) {
    FE += gamma[i][1] * log(pstar->pi[i]);
  }

  for (t=1; t<=T-1; t++) {
    for (i=1; i <= pstar->I; i++) {
      for (j=1; j <= pstar->I; j++) {
        FE += xi[i][j][t] * log(pstar->A[i][j]);
      }
        }
  }

  for (t=1; t<=T; t++) {
    for (i=1; i <= pstar->I; i++) {
      FE += gamma[i][t] * log(pstar->B[i][X[t]]);
        }
  }

/*  return -FE; */
  return FE;
}

/*           Classical_E_Step_Scaled(&pdhmm, T, X, alpha, beta, xi, gamma, scale);
 *
 * Purpose:  Perform the classical E-step with scaling.
 *
 *   Input:
 *  Output:  none
 *  Return:  none
 */

void Classical_E_step_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, 
              double **beta, double ***xi, double **gamma, double *scale)
{
  Forward_Scaled(pdhmm, T, X, alpha, scale);
  Backward_Scaled(pdhmm, T, X, beta, scale);

  Accumulate_xi_Scaled(pdhmm, T, X, alpha, beta, xi, scale);
  Normalize_xi(xi, pdhmm->I, pdhmm->M, T);
  
  Accumulate_gamma(pdhmm, T, X, alpha, beta, gamma);
  Normalize_gamma(gamma, pdhmm->I, T);
}


