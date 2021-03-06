/*
* H. Nahrstaedt - Aug 2010
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

/* IMPORTANT: uncomment the following line if you experience MATLAB chrashes */
/* instead of a normal error message when calling the function with bad syntax */
/* This bug is apparently restricted to GNU/Linux systems and to MATLAB versions prior to R2007a */
/*#define _ALT_MEXERRMSGTXT_*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "emd.h"



int c_stop_sifting(COMPLEX_T *, double *,extrema_t *,stop_t *,int,int);



/************************************************************************/
/*                                                                      */
/* MAIN FUNCTION                                                   */
/*                                                                      */
/************************************************************************/

//void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  int int_cemdc(char *fname)
{
  
  /* declarations */
  int i,n,nb_imfs,max_imfs,iteration_counter,stop_status,allocated_x,stop_EMD,nbphases;
  extrema_t ex;
  c_input_t input;
  envelope_t env;
  stop_t stop_params;
  double *x,*a;
  COMPLEX_T *y,*m,*z;
  c_imf_list_t list;
  FILE *fid;
      SciErr sciErr;
  CheckRhs (2,5);
  CheckLhs (1,2);
  /* get input data */
  input=c_get_input(pvApiCtx,Rhs,&sciErr);
  if(sciErr.iErr)
  {
    printError(&sciErr, 0);
    return 0;
  }
  n=input.n;
  nbphases=input.nbphases;
  max_imfs=input.max_imfs;
  stop_params=input.stop_params;
  allocated_x=input.allocated_x;
  x=input.x;
  y=input.y;
  
  /* initialisations */
  list=c_init_imf_list(n);
  z=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));
  m=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));
  a=(double *)malloc(n*sizeof(double));
  ex=c_init_extr(n+2*NBSYM);
  env=c_init_local_mean(n+2*NBSYM);
  
  /* MAIN LOOP */
  
  nb_imfs=0;
  stop_EMD=0;
  
  while ((!max_imfs || (nb_imfs < max_imfs)) && !stop_EMD) {
    
    /* initialisation */
    for (i=0;i<n;i++) z[i]=y[i];
    for (i=0;i<n;i++) m[i]=y[i];
    iteration_counter=0;
    
    stop_status = c_mean_and_amplitude(x,z,m,a,n,nbphases,&ex,&env);
   
    /* SIFTING LOOP */
    
    while (!stop_status && !c_stop_sifting(m,a,&ex,&stop_params,n,iteration_counter)) {
      
      /* subtract the local mean */
      #ifdef C99_OK
      for (i=0;i<n;i++) z[i]=z[i]-m[i];
      #else
      for (i=0;i<n;i++) {
        z[i].r=z[i].r-m[i].r;
        z[i].i=z[i].i-m[i].i;
      }
      #endif
      iteration_counter++;

      stop_status = c_mean_and_amplitude(x,z,m,a,n,nbphases,&ex,&env);
      
      
    }
    
    /* save current IMF into list if at least     */
    /* one sifting iteration has been performed */
    if (iteration_counter) {
      c_add_imf(&list,z,iteration_counter);
      nb_imfs++;
      #ifdef C99_OK
      for (i=0;i<n;i++) y[i]=y[i]-z[i];
      #else
      for (i=0;i<n;i++) {
        y[i].r=y[i].r-z[i].r;
        y[i].i=y[i].i-z[i].i;
      }
      #endif
    }
    else
      stop_EMD = 1;
    
  }
  
  /* save the residual into list */
  c_add_imf(&list,y,0);
  
  /* output into a MATLAB array */
  c_write_output(list,pvApiCtx,Rhs,Lhs);
  
   LhsVar(1) = Rhs + 1; 
   if (Lhs>1)
      LhsVar(2) = Rhs + 2; 
   PutLhsVar();
  /* free allocated memory */
  //if (allocated_x)
    free(x);
  free(y);
  free(m);
  free(a);
  c_free_local_mean(env);
  free(z);
  c_free_imf_list(list);
  c_free_extr(ex);
    return 0;
}

/************************************************************************/
/* STOP TEST FOR THE SIFTING LOOP                                    */
/************************************************************************/

int c_stop_sifting(COMPLEX_T *m, double *a,extrema_t *ex,stop_t *sp,int n, int counter) {
  int i,count;
  double tol,eps;
  tol = sp->tolerance*n;
  eps = sp->threshold;
  count = 0;
  if (counter >= MAX_ITERATIONS) return 1;
  for (i=0;i<n;i++) {
    if (CABS(m[i]) > eps*a[i]) if (++count>tol) return 0;
  }
  return 1;
}
