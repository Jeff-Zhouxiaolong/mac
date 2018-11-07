/*
* H. Nahrstaedt - Aug 2010
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#include "emd.h"

/************************************************************************/
/*                                                                      */
/* GET INPUT DATA                                                       */
/*                                                                      */
/************************************************************************/

c_input_t c_get_input(void* _pvCtx,int nrhs, SciErr* sciErr) {
  c_input_t input;
  int n,i,nY;
  double *x,*x_temp,*ry_temp,*iy_temp,*third,*fourth,*fifth;
  COMPLEX_T *y;
  
  int* piAddressVar     = NULL;
  int iType = 0;  
  int mVar = 0, nVar = 0;
  
  input.stop_params.threshold = DEFAULT_THRESHOLD;
  input.stop_params.tolerance = DEFAULT_TOLERANCE;
  input.allocated_x=0;

  input.max_imfs=0;
  input.nbphases=DEFAULT_NBPHASES;
  
  /* argument checking*/
//   if (nrhs>5)
//     mexErrMsgTxt("Too many arguments");
//   if (nrhs<2)
//     mexErrMsgTxt("Not enough arguments");
//   if (nlhs>2)
//     mexErrMsgTxt("Too many output arguments");

  /* get Address of inputs */
  *sciErr = getVarAddressFromPosition(_pvCtx, 1, &piAddressVar);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check input type */
  *sciErr = getVarType(_pvCtx, piAddressVar, &iType);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if ( iType != sci_matrix )
  {
    Scierror(999," T must be either empty or a double precision real vector.\n");
    sciErr->iErr=999;
    return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&x_temp);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check size */
  if (( (mVar !=1) && (nVar !=1 ) ) && !(mVar==0 && nVar==0))
  {
  	Scierror(999,"T must be either empty or a double precision real vector.\n");
	sciErr->iErr=999;
  	return input;
  }
     nY=GREATER(mVar,nVar);


//   if (!mxIsEmpty(prhs[0]))
//     if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
//     mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) ||
//     (mxGetNumberOfDimensions(prhs[0]) > 2))
//       mexErrMsgTxt("X must be either empty or a double precision real vector.");
     
       /* get Address of inputs */
  *sciErr = getVarAddressFromPosition(_pvCtx, 2, &piAddressVar);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check input type */
  *sciErr = getVarType(_pvCtx, piAddressVar, &iType);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if (( iType != sci_matrix ) || ! isVarComplex(_pvCtx, piAddressVar))
  {
    Scierror(999,"X must be a double precision complex vector.\n");
    sciErr->iErr=999;
    return input;
  }
     
     
   *sciErr = getComplexMatrixOfDouble(_pvCtx, piAddressVar, &mVar, &nVar, &ry_temp, &iy_temp);

  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if ( (mVar !=1) && (nVar !=1 ) ) 
  {
  	Scierror(999,"X must be  a double precision complex vector..\n");
	sciErr->iErr=999;
  	return input;
  }
    n=GREATER(mVar,nVar);  
     
    
  
//   if (!mxIsNumeric(prhs[1]) || !mxIsComplex(prhs[1]) ||
//   mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) ||/* length of vector x */
//   (mxGetNumberOfDimensions(prhs[1]) > 2))
//     mexErrMsgTxt("Y must be a double precision complex vector.");
  
  /* input reading: x and y */
//   n=GREATER(mxGetN(prhs[1]),mxGetM(prhs[1])); /* length of vector x */
//   if (mxIsEmpty(prhs[0])) {
  x = (double *)malloc(n*sizeof(double));
if (nY==0){ 
    input.allocated_x = 1;
    for(i=0;i<n;i++) x[i] = i;
    nY=n;
  }
  else{ 
   for (i=0;i<n;i++)
      x[i]=x_temp[i];
  }
    
    
//     x=mxGetPr(prhs[0]);
//   ry_temp=mxGetPr(prhs[1]);
//   iy_temp=mxGetPi(prhs[1]);

   /* check size */
  if ( (n !=nY) ) 
  {
  	Scierror(999,"X and T must have the same length.\n");
	sciErr->iErr=999;
  	return input;
  }
  
  /* third argument */
  if (nrhs>=3) {
    
      /* get Address of inputs */
  *sciErr = getVarAddressFromPosition(_pvCtx, 3, &piAddressVar);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check input type */
  *sciErr = getVarType(_pvCtx, piAddressVar, &iType);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if ( iType != sci_matrix )
  {
    Scierror(999,"STOP must be a real vector of 1 or 2 elements.\n");
    sciErr->iErr=999;
    return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&third);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
  i = GREATER(mVar,nVar);
    /* check size */
  if ( i>2 ) 
  {
  	Scierror(999,"STOP must be a real vector of 1 or 2 elements.\n");
	sciErr->iErr=999;
  	return input;
  }
    

  
//     if(!mxIsEmpty(prhs[2])) {
//       if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])
//       || !mxIsDouble(prhs[2]) || (mxGetN(prhs[2])!=1 && mxGetM(prhs[2])!=1))
//         mexPrintf("STOP must be a real vector of 1 or 2 elements");
//       i = GREATER(mxGetN(prhs[2]),mxGetM(prhs[2]));
//       if (i>2)
//         mexErrMsgTxt("STOP must be a vector of 1 or 2 elements");
//       third=mxGetPr(prhs[2]);
      
    if (i>0) {    
      switch (i) {
        case 1 : {
          input.stop_params.threshold=third[0];
          break;
        }
        case 2 : {
          input.stop_params.threshold=third[0];
          input.stop_params.tolerance=third[1];
        }
      }
      /* input checking */
      if (input.stop_params.threshold <= 0)  {
        Scierror(999,"threshold must be a positive number");
        sciErr->iErr=999;
	return input;
    }
      if (input.stop_params.threshold >= 1){
         Scierror(999,"threshold should be lower than 1");
	 sciErr->iErr=999;
	 return input;
      }
      if (input.stop_params.tolerance < 0 || input.stop_params.tolerance >= 1){
         Scierror(999,"tolerance must be a real number in [O,1]");
	 sciErr->iErr=999;
	 return input;
      }
    }
  }
  
  
  /* fourth argument */
  if (nrhs>=4) {
    
    
    /* get Address of inputs */
  *sciErr = getVarAddressFromPosition(_pvCtx, 4, &piAddressVar);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check input type */
  *sciErr = getVarType(_pvCtx, piAddressVar, &iType);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if ( iType != sci_matrix )
  {
    Scierror(999,"NB_IMFS must be a positive integer.\n");
    sciErr->iErr=999;
	 return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&fourth);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
      if ( ( (mVar !=1) || (nVar !=1) )&& (mVar>0 && nVar>0))
  {
        Scierror(999,"NB_IMFS must be a positive integer.\n");
       sciErr->iErr=999;
	 return input;
  }
    
  if (mVar>0 && nVar>0){
     if ((unsigned int)fourth[0] != fourth[0]){
        Scierror(999,"NB_IMFS must be a positive integer");
	sciErr->iErr=999;
	 return input;
     }
      input.max_imfs=(int)fourth[0];
   
    
//     if (!mxIsEmpty(prhs[3])) { /* if empty -> do nothing */
//       if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3])
//       || !mxIsDouble(prhs[3]) || mxGetN(prhs[3])!=1 || mxGetM(prhs[3])!=1)
//         mexErrMsgTxt("NB_IMFS must be a positive integer");
//       fourth=*mxGetPr(prhs[3]);
//       if ((unsigned int)fourth != fourth)
//         mexErrMsgTxt("NB_IMFS must be a positive integer");
//       input.max_imfs=(int)fourth;
 
    }
  }
  
  /* fifth argument */
  if (nrhs==5) {
    
     /* get Address of inputs */
  *sciErr = getVarAddressFromPosition(_pvCtx, 5, &piAddressVar);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    /* check input type */
  *sciErr = getVarType(_pvCtx, piAddressVar, &iType);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
     if ( iType != sci_matrix )
  {
    Scierror(999,"NBPHASES must be a positive integer.\n");
    sciErr->iErr=999;
	 return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&fifth);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
      if ( ( (mVar !=1) || (nVar !=1) )&& (mVar>0 && nVar>0))
  {
        Scierror(999,"NBPHASES must be a positive integer.\n");
       sciErr->iErr=999;
	 return input;
  }
    if (mVar>0 && nVar>0){
     if ((int)fifth[0] != fifth[0]){
        Scierror(999,"NBPHASES must be a positive integer");
	sciErr->iErr=999;
	 return input;
     }
      input.nbphases = (int)fifth[0];

//     if(!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4])
//     || mxGetN(prhs[4])!=1 || mxGetM(prhs[4])!=1)
//       mexErrMsgTxt("NBPHASES must be a positive integer");
//     fifth=*mxGetPr(prhs[4]);
//     if ((int)fifth != fifth)
//       mexErrMsgTxt("NBPHASES must be a positive integer");
//     input.nbphases = (int)fifth;
      }
  }
  
  /* more input checking */
//   if (!input.allocated_x && SMALLER(mxGetN(prhs[0]),mxGetM(prhs[0]))!=1 ||
//   SMALLER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=1)
//     mexErrMsgTxt("X and Y must be vectors");
//   if (GREATER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=n)
//     mexErrMsgTxt("X and Y must have the same length");
  i=1;
  while (i<n && x[i]>x[i-1]) i++;
  if (i<n) {
    Scierror(999,"Values in T must be non decreasing");
    sciErr->iErr=999;
	 return input;
  }
  
  /* copy vector y to avoid erasing input data */
  y=(COMPLEX_T *)malloc(n*sizeof(COMPLEX_T));
  #ifdef C99_OK
  for (i=0;i<n;i++) y[i]=ry_temp[i]+I*iy_temp[i];
  #else
  for (i=0;i<n;i++) {
    y[i].r=ry_temp[i];
    y[i].i=iy_temp[i];
  }
  #endif
  
  input.n=n;
  input.x=x;
  input.y=y;
  return input;
}


/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF THE LIST                                           */
/*                                                                      */
/************************************************************************/

c_imf_list_t c_init_imf_list(int n) {
  c_imf_list_t list;
  list.first=NULL;
  list.last=NULL;
  list.n=n;
  list.m=0;
  return list;
}


/************************************************************************/
/*                                                                      */
/* ADD AN IMF TO THE LIST                                               */
/*                                                                      */
/************************************************************************/

void c_add_imf(c_imf_list_t *list,COMPLEX_T *p,int nb_it) {
  COMPLEX_T *v=(COMPLEX_T *)malloc(list->n*sizeof(COMPLEX_T));
  int i;
  c_imf_t *mode=(c_imf_t *)malloc(sizeof(c_imf_t));
  for (i=0;i<list->n;i++) v[i]=p[i];
  mode->pointer=v;
  mode->nb_iterations=nb_it;
  mode->next=NULL;
  if (!list->first) {
    list->first=mode;
  } else {
    (list->last)->next=mode;
  }
  list->last=mode;
  list->m++;
}


/************************************************************************/
/*                                                                      */
/* FREE MEMORY ALLOCATED FOR THE LIST                                   */
/*                                                                      */
/************************************************************************/

void c_free_imf_list(c_imf_list_t list) {
  c_imf_t *current=list.first, *previous;
  while (current) {
    previous=current;
    current=current->next;
    free(previous->pointer);
    free(previous);
  }
}


/************************************************************************/
/*                                                                      */
/* OUTPUT INTO MATLAB ARRAY                                             */
/*                                                                      */
/************************************************************************/

void c_write_output(c_imf_list_t list,void* _pvCtx,int rhs,int lhs) {
  double *rout,*iout,*out2;
  c_imf_t *current;
  int i=0,j,m=list.m,n=list.n;
  rout =  (double*)malloc(sizeof(double) * m*n);
  iout =  (double*)malloc(sizeof(double) * m*n);
  if (lhs>1)
    out2 =  (double*)malloc(sizeof(double) * 1*(m-1));
//   plhs[0]=mxCreateDoubleMatrix(m,n,mxCOMPLEX);
//   rout=mxGetPr(plhs[0]);
//   iout=mxGetPi(plhs[0]);
//   plhs[1]=mxCreateDoubleMatrix(1,m-1,mxCOMPLEX);
//   out2=mxGetPr(plhs[1]);

  for (current=list.first;current;current=current->next) {
    for (j=0;j<n;j++) {
      *(rout+j*m+i)=CREAL(current->pointer[j]);
      *(iout+j*m+i)=CIMAG(current->pointer[j]);
    }
    if ((i<m-1)&& lhs>1)
          *(out2+i)=current->nb_iterations;
    i++;
  }
  
  createComplexMatrixOfDouble(_pvCtx, rhs + 1,m,n,rout,iout);
     free(rout);   free(iout);
   if (lhs>1){  
    createMatrixOfDouble(_pvCtx, rhs + 2, 1, m-1, out2);
    free(out2);
  }
  
  
}
