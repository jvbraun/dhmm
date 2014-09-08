/*      Program:  algorithms.c
 *      Version:  0.1
 *         Date:  1 February 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Algorithms to work with hidden Markov models.
 */

static char rcsid[] = "$Id: algorithms.nw,v 1.4 1999/07/10 09:03:08 Jerome Braun Exp $";

/*********************************************************************
 * INCLUDES 
 ********************************************************************/

#include <stdlib.h>             /* we want the rand function */
#include <math.h>               /* we want the log function */
#include <stdio.h>              /* we want to see what is going on... */
#include "dhmm.h"
#include "utilities.h"

/*********************************************************************
 * DEFINES 
 ********************************************************************/

#define NCYCLE 1000              /* number of cycles to iterate */
#define CRITERION 1e-03         /* ratio convergence criterion */
#define DELTA_OLD 1e-70         /* some starting number */

/*********************************************************************
 * ALGORITHMS 
 ********************************************************************/

/*           Baum_Welch(dhmm, T, X, alpha, beta, gamma, xi);
 *
 * Purpose:  Carry out the classical Baum-Welch algorithm on one sequence
 *           and put the resulting estimate into dhmm.
 *
 *  Method:  The classical Baum-Welch algorithm is an instance of the E-M
 *           algorithm, iterative expectation-maximization of the maximum
 *           likelihood.  This version implements scaling.
 *
 *   Input:  pointer to initialized DHMM structure for model,
 *           integer data array length, pointer to data array with data, 
 *           pointers to arrays for alpha, beta, gamma, and xi
 *  Output:  none
 *  Return:  none
 */

void Baum_Welch(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi)
{
  int cycle;                      /* cycle iterator */

  double *scale;                  /* scaling array */
  double LL, LL_old;              /* log-likelihood tracking */
  double delta, delta_old;        /* convergence criteria tracking */

  scale = array_double_1d(1, T);  /* initialize scaling array */
  delta_old = DELTA_OLD;          /* initialize the change */

  /* Initial M-step: handled already by initializing the hidden Markov structure */

  /* Initial E-step */
  Classical_E_step_Scaled(pdhmm, T, X, alpha, beta, xi, gamma, scale);

  LL = Log_Likelihood_Scaled(scale, T);
  LL_old = LL;  /* keep the log-likelihood for comparative purposes */

  /* now carry out the iterative E-M algorithm, a.k.a. Baum-Welch or */
  /* free energy maximization                                        */

  for (cycle=1; cycle<=NCYCLE; cycle++) {

    /* M-step */
    Classical_M_step(pdhmm, 1e-06, T, X, gamma, xi);

    /* E-step */
    Classical_E_step_Scaled(pdhmm, T, X, alpha, beta, xi, gamma, scale);

    /* use the log-likelihood to measure convergence */
    LL = Log_Likelihood_Scaled(scale, T);
    delta = LL - LL_old;
    
    if (delta / delta_old <= CRITERION) break; /* converged */

    /* not converged, keep going */
    LL_old = LL;
    delta_old = delta;

  } /* cycle again */

  /* clean up */
  free_array_double_1d(scale, 1, T);
}

/*           Penalized_Baum_Welch(dhmm, prior, T, X, alpha, beta, gamma, xi);
 *
 * Purpose:  Carry out the penalized Baum-Welch algorithm on one sequence
 *           and put the resulting estimate into dhmm.
 *
 *  Method:  The penalized version of the Baum_Welch algorithm corresponds to
 *           maximizing the posterior distribution given the prior.  
 *           This version implements scaling.
 *
 *   Input:  pointers to initialized DHMM structures for model and prior,
 *           integer data array length, pointer to data array with data, 
 *           pointers to arrays for alpha, beta, gamma, and xi
 *  Output:  none
 *  Return:  none
 */

void Penalized_Baum_Welch(DHMM *pdhmm, DHMM *pprior, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi)
{
  int cycle;                      /* cycle iterator */

  double *scale;                  /* scaling array */
  double LP, LP_old;              /* log-posterior-probability tracking */
  double delta, delta_old;        /* convergence criteria tracking */

  scale = array_double_1d(1, T);  /* initialize scaling array */
  delta_old = DELTA_OLD;          /* initialize the change */

  /* Initial M-step: handled already by initializing the hidden Markov structure */

  /* Initial E-step */
  Classical_E_step_Scaled(pdhmm, T, X, alpha, beta, xi, gamma, scale);

  /* in this case, the same computation that gave us the log-likelihood
   * now gives us the log-posterior-probability
   */
  LP = Log_Likelihood_Scaled(scale, T);
  LP_old = LP;  /* keep for comparative purposes */

  /* now carry out the iterative E-M algorithm, a.k.a. Baum-Welch or */
  /* free energy maximization                                        */

  for (cycle=1; cycle<=NCYCLE; cycle++) {

    /* M-step */
    Penalized_M_step(pdhmm, pprior, T, X, gamma, xi);

    /* E-step */
    Classical_E_step_Scaled(pdhmm, T, X, alpha, beta, xi, gamma, scale);

    /* use the posterior direct probability to measure convergence */
    LP = Log_Likelihood_Scaled(scale, T);
    delta = LP - LP_old;

    if (delta / delta_old <= CRITERION) break;  /* converged */
    
    /* not converged, keep going */
    LP_old = LP;
    delta_old = delta;

  } /* cycle again */

  /* clean up */
  free_array_double_1d(scale, 1, T);
}

/*           Approximate_Bayes(&posterior, &prior, T, X, alpha, beta, gamma, xi);
 *
 * Purpose:  Carry out the "ensemble learning algorithm" of David J.C. MacKay
 *           on one sequence and put the resulting approximation of the
 *           posterior distribution into posterior.
 *
 *  Method:  The "ensemble learning algorithm" provides an elegant 
 *           approximation to the full posterior distribution of the hidden
 *           Markov model when the prior is taken to be a product of
 *           independent Dirichlet distributions.
 * 
 *           This version implements scaling.
 *
 *   To do:  1.  Determine whether free energy is the best criterion.
 *
 *   Input:  pointers to initialized DHMM structures for model and prior,
 *           integer data array length, pointer to data array with data, 
 *           pointers to arrays for alpha, beta, gamma, and xi
 *  Output:  none
 *  Return:  none
 */

void Approximate_Bayes(DHMM *pposterior, DHMM *pprior, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi)
{
  int cycle;                /* cycle iterator */

  double *scale;            /* scaling array */
  double FE, FE_old;        /* free energy tracking */
  double delta, delta_old;  /* convergence criterion tracking */

  DHMM star;                /* DHMM structure for A^*, B^*, pi^* */

  new_DHMM(&star, pprior->I, pprior->M);   /* allocate the star matrix */
  scale = array_double_1d(1, T);           /* initialize scaling array */

  delta_old = DELTA_OLD;                   /* initialize the change */

  /* Initial M-step: handled already by initializing the estimate */
  Calculate_star(pposterior, &star);

  /* Initial E-step */
  Classical_E_step_Scaled(&star, T, X, alpha, beta, xi, gamma, scale);


  FE = Free_Energy(&star, pprior, X, T, gamma, xi); 
  FE_old = FE;  /* keep the free energy for comparative purposes */

  /* now carry out the iterative E-M algorithm, a.k.a. */
  /* free energy maximization                                        */

  for (cycle=1; cycle<=NCYCLE; cycle++) {

    /* M-step */
    Approximate_Bayes_M_step(pposterior, pprior, T, X, gamma, xi);
    Calculate_star(pposterior, &star);

    /* E-step */
    Classical_E_step_Scaled(&star, T, X, alpha, beta, xi, gamma, scale);

    FE = Free_Energy(&star, pprior, X, T, gamma, xi); 
   /* use the free energy as a measure of convergence */
   delta = FE_old - FE;

   /* now we need to work through two cycles before comparison */
   if (cycle>2 && delta <= CRITERION) break; /* converged */
   
   /* not converged, keep going */
    FE_old = FE;
    delta_old = delta;

  } /* cycle */

  /* clean up */
  free_array_double_1d(scale, 1, T);
  free_DHMM(&star);
}  

