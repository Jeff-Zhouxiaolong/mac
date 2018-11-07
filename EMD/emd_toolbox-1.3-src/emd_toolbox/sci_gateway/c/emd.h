
/*
* H. Nahrstaedt - Aug 2010
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/



#ifndef EMD_H
#define EMD_H 1

// #include <stack-c.h>
#include "api_scilab.h"
#include "Scierror.h"
#include "MALLOC.h"

#ifdef C99_OK
  #include <complex.h>
  #define COMPLEX_T double complex
  #define CREAL creal
  #define CIMAG cimag
  #define CABS cabs
#else
  #define COMPLEX_T complex_data_t
  #define CREAL emd_creal
  #define CIMAG emd_cimag
  #define CABS emd_cabs
#endif

#define GREATER(A,B) ((A) >= (B) ? (A) : (B))
#define SMALLER(A,B) ((A) <  (B) ? (A) : (B))

#define SQUARE(A) (A*A)
#define CUBE(A) (A*A*A)

#define STOP_DEFAULT {.threshold = 0.05, .tolerance = 0.05}
#define DEFAULT_THRESHOLD 0.05
#define DEFAULT_TOLERANCE 0.05
#define MAX_ITERATIONS 1000
#define LIM_GMP 30000
#define NBSYM 2
#define DEFAULT_NBPHASES 4

//emdc_fix
#define DEFAULT_NB_ITERATIONS 10

/* structure used to store envelopes and temporary data */
// typedef struct {
//   int n;
//   double *e_min;
//   double *e_max;
//   double *tmp1;
//   double *tmp2;
// } envelope_t;


typedef struct {
    double threshold;
    double tolerance;
} stop_t;

/* structure used to store input data */
typedef struct {
  int n;
  int max_imfs;
  int allocated_x;
  int nb_iterations;
  double *x;
  double *y;
  stop_t stop_params;
} input_t;


/* structure used to store an IMF and the associated number of iterations */
typedef struct i {
  int nb_iterations;
  double *pointer;
  struct i *next;
} imf_t;

/* structure of the IMF list */
typedef struct {
  imf_t *first;
  imf_t *last;
  int m; 
  int n; 
} imf_list_t;

typedef struct {
  double r;
  double i;
} complex_data_t;

/* structure used to store input data */
typedef struct {
  int n;
  int max_imfs;
  int nb_iterations;
  int allocated_x;
  int nbphases;
  double *x;
  COMPLEX_T *y;
  stop_t stop_params;
} c_input_t;


/* structure used to store an IMF and the associated number of iterations */
typedef struct c_i {
  int nb_iterations;
  COMPLEX_T *pointer;
  struct c_i *next;
} c_imf_t;

/* structure of the IMF list */
typedef struct {
  c_imf_t *first;
  c_imf_t *last;
  int m; 
  int n; 
} c_imf_list_t;




/* structure used to store envelopes and temporary data */
typedef struct {
  int n;
  double *re_min;
  double *ie_min;
  double *re_max;
  double *ie_max;
  double *e_min;
  double *e_max; 
  double *tmp1;
  double *tmp2;
} envelope_t;

// typedef struct {
//   int n_min;
//   int n_max;
//   double *x_min;
//   double *y_min;
//   double *x_max;
//   double *y_max;
// } extrema_t;
/* structure used to store the local extrema of the signal */
typedef struct {
  int n_min;
  int n_max;
  int *ind_min;
  int *ind_max;
  double *x_min;
  double *ry_min;
  double *iy_min;
  double *y_min;  
  double *x_max;
  double *ry_max;
  double *iy_max;
   double *y_max; 
} extrema_t;

extrema_t init_extr(int);
void extr(double *,double *,int,extrema_t *);
void free_extr(extrema_t);
void boundary_conditions(double *,double *,int,extrema_t *);

void interpolation(double *,double *,double *,int,double *,int,double *, double *);

envelope_t init_local_mean(int);
void free_local_mean(envelope_t);
int mean(double *,double *,double *,int,extrema_t *,envelope_t *);
int mean_and_amplitude(double *,double *,double *,double *,int,extrema_t *,envelope_t *);


input_t get_input(void* _pvCtx,int nrhs, SciErr* sciErr);
input_t get_input_fix(void* _pvCtx,int nrhs, SciErr* sciErr);
imf_list_t init_imf_list(int);
void add_imf(imf_list_t *,double *,int);
void free_imf_list(imf_list_t);
void write_output(imf_list_t list,void* _pvCtx,int rhs,int lhs);

//complex versions
#ifndef C99_OK
double CREAL(COMPLEX_T);
double CIMAG(COMPLEX_T);
double CABS(COMPLEX_T);
double crealeiphi(double,COMPLEX_T);
#endif


extrema_t c_init_extr(int);
void c_extr(double *,COMPLEX_T *,double,int,extrema_t *);
void c_free_extr(extrema_t);
void c_boundary_conditions(double *,COMPLEX_T *,double,int,extrema_t *);

envelope_t c_init_local_mean(int);
int c_mean_and_amplitude(double *,COMPLEX_T *,COMPLEX_T *,double *,int,int,extrema_t *,envelope_t *);
int c_mean(double *,COMPLEX_T *,COMPLEX_T *,int,int,extrema_t *,envelope_t *);
void c_free_local_mean(envelope_t);

envelope_t c2_init_local_mean(int);
int c2_mean_and_amplitude(double *,COMPLEX_T *,COMPLEX_T *,double *,int,int,extrema_t *,envelope_t *);
int c2_mean(double *,COMPLEX_T *,COMPLEX_T *,int,int,extrema_t *,envelope_t *);
void c2_free_local_mean(envelope_t);

c_input_t c_get_input(void* _pvCtx,int nrhs, SciErr* sciErr);
c_input_t c_get_input_fix(void* _pvCtx,int nrhs, SciErr* sciErr);
c_imf_list_t c_init_imf_list(int);
void c_add_imf(c_imf_list_t *,COMPLEX_T *,int);
void c_free_imf_list(c_imf_list_t);
void c_write_output(c_imf_list_t,void* _pvCtx,int rhs,int lhs);

#endif /* EMD_H */