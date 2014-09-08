/*      Program:  dhmm.h
 *      Version:  0.1
 *         Date:  1 February 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Typedefs and function prototypes for working
 *                with discrete hidden Markov models.
 */

/*********************************************************************
 * TYPEDEFS 
 ********************************************************************/

/* Structure to contain the parameters of a discrete hidden Markov
 * model.  (Nice structure due to Tapas Kanungo.) 
 */

typedef struct {
  int I;                /* number of states */
  int M;                /* number of observations */
  double **A;   /* state transition matrix */
  double **B;   /* observation emission matrix */
  double *pi;   /* initial state distribution */
} DHMM;

/*********************************************************************
 * FUNCTION PROTOTYPES 
 ********************************************************************/

/* utilities for working with DHMM structures */
void initialize_DHMM(unsigned long *seed, DHMM *pdhmm);
void free_DHMM(DHMM *pdhmm);
void read_DHMM(FILE *fp, DHMM *pdhmm);
void write_DHMM(FILE *fp, DHMM *pdhmm);
void DHMM_deviate(unsigned long *seed, DHMM *pprior, DHMM *pdhmm);
void copy_DHMM(DHMM *poriginal, DHMM *pcopy);
void new_DHMM(DHMM *pdhmm, int I, int M);

/* utilities for working with sequences to be analyzed */
void read_sequence(FILE *fp, int *pT, int **pX);
void write_sequence(FILE *fp, int T, int *X);
void generate_sequence(unsigned long *seed, DHMM *pdhmm, int T, int **pX);

/* subroutines to carry out various steps of the algorithms */
void Forward(DHMM *pdhmm, int T, int *X, double **alpha);
void Forward_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, double *scale);
void Backward(DHMM *pdhmm, int T, int *X, double **beta);
void Backward_Scaled(DHMM *pdhmm, int T, int *X, double **beta, double *scale);

void Normalize_xi(double ***xi, int I, int M, int T);
void Normalize_gamma(double **gamma, int I, int T);
void Accumulate_xi(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double ***xi);
void Accumulate_xi_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double ***xi, double *scale);
void Accumulate_gamma(DHMM *pposterior, int T, int *X, double **alpha, double **beta, double **gamma);
void Calculate_star(DHMM *pposterior, DHMM *pstar);
void Classical_E_step_Scaled(DHMM *pdhmm, int T, int *X, double **alpha, 
              double **beta, double ***xi, double **gamma, double *scale);

double Free_Energy(DHMM *pdhmm, DHMM *pprior, int *X, int T, double **gamma, double ***xi);
double Log_Likelihood(double **alpha, int I, int T);
double Log_Likelihood_Scaled(double *scale, int T);

void Classical_M_step (DHMM *pdhmm, double epsilon, int T, int *X, double **gamma, double ***xi);
void Penalized_M_step (DHMM *pdhmm, DHMM *pprior, int T, int *X, double **gamma, double ***xi);
void Approximate_Bayes_M_step (DHMM *pposterior, DHMM *pprior, int T, int *X, double **gamma, double ***xi);

/* complete algorithms */
void Baum_Welch(DHMM *pdhmm, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi);
void Penalized_Baum_Welch(DHMM *pdhmm, DHMM *pprior, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi);
void Approximate_Bayes(DHMM *pposterior, DHMM *pprior, int T, int *X, double **alpha, double **beta, double **gamma, double ***xi);
