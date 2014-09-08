/*      Program:  utilities.h
 *      Version:  0.1
 *         Date:  26 February 1999
 *       Author:  Jerome V. Braun
 *       E-mail:  jerome.braun@kmri.com
 * Organization:  Kings Mountain Research
 *              
 *      Purpose:  Function prototypes for the utility functions
 *                in utilities.h 
 */

/* Utility function prototypes */

void utilities_error();

/* Array memory allocation function prototypes */

int *array_int_1d();

double *array_double_1d();
double **array_double_2d();
double ***array_double_3d();

void free_array_int_1d();
void free_array_double_1d();
void free_array_double_2d();
void free_array_double_3d();

/* Mathematical function prototypes */

double digamma();

/* Random number generator wrapper functions used in an
 * effort to make changes at least a bit easier if new
 * random number generators are desired
 */

void set_seed(unsigned long *seed, unsigned long k);
void save_seed(unsigned long *seed);
double uniform_deviate(unsigned long *seed);
int multinomial_deviate(unsigned long *seed, double *p, int first, int last);
double gamma_deviate(unsigned long *seed, double a);
void Dirichlet_deviate(unsigned long *seed, double *u, int i1, int i2, double *p); 

/* Win_rand 1.0 function prototypes which are wrapped
 * by the above functions
 */

double drand(unsigned long *seed);
double gll(unsigned long *seed, double a);
