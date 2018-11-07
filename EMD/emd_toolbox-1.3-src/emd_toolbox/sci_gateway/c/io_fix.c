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

input_t get_input_fix(void* _pvCtx,int nrhs, SciErr* sciErr){
  input_t input;
  int n,i,nY;
  double *x,*x_temp,*y,*y_temp,*third,*fourth;
  
  int* piAddressVar     = NULL;
  int iType = 0;  
  int mVar = 0, nVar = 0;
  
  input.nb_iterations=DEFAULT_NB_ITERATIONS;
  input.allocated_x=0;

  input.max_imfs=0;

  
//   /* argument checking*/
//   if (nrhs>4)
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
    Scierror(999,"T must be either empty or a double precision real vector.\n");
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
     if ( iType != sci_matrix )
  {
    Scierror(999,"X must be  a double precision real vector.\n");
    sciErr->iErr=999;
    return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&y_temp);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
  if ( (mVar !=1) && (nVar !=1 ) ) 
  {
  	Scierror(999,"X must be  a double precision real vector..\n");
	sciErr->iErr=999;
  	return input;
  }
  
  
    n=GREATER(mVar,nVar);
    
    
    
  
//   if (!mxIsEmpty(prhs[0]))
//     if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
//     mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) ||
//     (mxGetNumberOfDimensions(prhs[0]) > 2))
//       mexErrMsgTxt("X must be either empty or a double precision real vector.");
//   
//   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
//   mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1])  ||
//   (mxGetNumberOfDimensions(prhs[1]) > 2))
//     mexErrMsgTxt("Y must be a double precision real vector.");
  
  /* input reading: x and y */
  //n=GREATER(mxGetN(prhs[1]),mxGetM(prhs[1])); /* length of vector x */
  //if (mxIsEmpty(prhs[0])) {
   x = (double *)malloc(n*sizeof(double));
 if (nY==0)  {
    input.allocated_x = 1;
    for(i=0;i<n;i++) x[i] = i;
    nY=n;
  }
  else
     for (i=0;i<n;i++) x[i]=x_temp[i];
//     x=mxGetPr(prhs[0]);
//   y_temp=mxGetPr(prhs[1]);
     
     
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
    Scierror(999,"The number of iterations must be a positive integer.\n");
    sciErr->iErr=999;
    return input;
  }
    *sciErr = getMatrixOfDouble(_pvCtx, piAddressVar,&mVar,&nVar,&third);
  if(sciErr->iErr)
  {
    printError(sciErr, 0);
    return input;
  }
    
     if ((  (mVar !=1) || (nVar !=1) )&& (mVar>0 && nVar>0))
  {
        Scierror(999,"The number of iterations must be a positive integer.\n");
       sciErr->iErr=999;
	 return input;
  }
    
  if (mVar>0 && nVar>0){
     if ((int)third[0] != third[0]){
        Scierror(999,"The number of iterations must be a positive integer");
	sciErr->iErr=999;
	 return input;
     }
       input.nb_iterations = (int)third[0];
  }
  }
    
    
//     if (!mxIsEmpty(prhs[2])) { /* if empty -> do nothing */
//       if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])
//       || !mxIsDouble(prhs[2]) || mxGetN(prhs[2])!=1 || mxGetM(prhs[2])!=1)
//         mexErrMsgTxt("The number of iterations must be a positive integer");
//       third=*mxGetPr(prhs[2]);
//       if((int)third != third)
//         mexErrMsgTxt("The number of iterations must be a positive integer");
//       input.nb_iterations = (int)third;
//     }
  
  
  
    /* fourth argument */
  if (nrhs==4) {
    
    
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
  }
    
    
    
    
    
    
//     if (!mxIsEmpty(prhs[3])) { /* if empty -> do nothing */
//       if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3])
//       || !mxIsDouble(prhs[3]) || mxGetN(prhs[3])!=1 || mxGetM(prhs[3])!=1)
//         mexErrMsgTxt("NB_IMFS must be a positive integer");
//       fourth=*mxGetPr(prhs[3]);
//       if ((unsigned int)fourth != fourth)
//         mexErrMsgTxt("NB_IMFS must be a positive integer");
//       input.max_imfs=(int)fourth;
//     }
  }
  
  /* more input checking */
//   if (!input.allocated_x && SMALLER(mxGetN(prhs[0]),mxGetM(prhs[0]))!=1 ||
//   SMALLER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=1)
//     mexErrMsgTxt("X and Y must be vectors");
//   if (GREATER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=n)
//     mexErrMsgTxt("X and Y must have the same length");
  i=1;
  while (i<n && x[i]>x[i-1]) i++;
  if (i<n){
     Scierror(999,"Values in T must be non decreasing");
    sciErr->iErr=999;
	 return input;
  }
  
  /* copy vector y to avoid erasing input data */
  y=(double *)malloc(n*sizeof(double));
  for (i=0;i<n;i++) y[i]=y_temp[i];
  
  input.n=n;
  input.x=x;
  input.y=y;
  return input;
}





/************************************************************************/
/*                                                                      */
/* OUTPUT INTO MATLAB ARRAY                                             */
/*                                                                      */
/************************************************************************/

// void write_output(imf_list_t list,mxArray *plhs[]) {
//   double *out1,*out2;
//   imf_t *current;
//   int i=0,j,m=list.m,n=list.n;
//   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
//   out1=mxGetPr(plhs[0]);
//   plhs[1]=mxCreateDoubleMatrix(1,m-1,mxREAL);
//   out2=mxGetPr(plhs[1]);
//   for (current=list.first;current;current=current->next) {
//     for (j=0;j<n;j++) *(out1+j*m+i)=current->pointer[j];
//     if (i<m-1) *(out2+i)=current->nb_iterations;
//     i++;
//   }
// }
