function [IMF,NB_ITERATIONS]=cemdc_fix(T,X,NB_ITERATIONS,MAX_IMFS,NDIRS);
//  bivariate Empirical Mode Decomposition, first algorithm
//  Calling Sequence
// [IMF,NB_ITERATIONS]=cemdc_fix(T,X);
// [IMF,NB_ITERATIONS]=cemdc_fix([],X);
// [IMF,NB_ITERATIONS]=cemdc_fix(T,X,NB_ITERATIONS);
// [IMF,NB_ITERATIONS]=cemdc_fix(T,X,NB_ITERATIONS,MAX_IMFS);
// [IMF,NB_ITERATIONS]=cemdc_fix(T,X,NB_ITERATIONS,MAX_IMFS,NDIRS);
// Parameters
// inputs:	
//       - T: sampling times. If T=[], the signal is assumed uniformly sampled.1xN time instants
//       - X: analyzed signal, 1xN signal data
//       - NB_ITERATIONS: number of sifting iterations to be performed to   extract each IMF. If NB_ITERATIONS is empty or unspecified, 10 iterations   are performed by default.  Note: The effective number of sifting iterations might be less  than NB_ITERATIONS for the last modes if the sifting process has to be stopped because of a lack of extrema.
//       - MAX_IMFS: maximum number of IMFs to be extracted. If MAX_IMFS is  zero, empty or unspecified, the default behavior is to extract as   many IMFs as possible.
//       - NDIRS: number of directions used to compute the local mean.  If unspecified, the default value is 4. TODO: the actual number of directions (according to [1]) is 2*NDIRS
// outputs: 
//              - IMF: intrinsic mode functions (IMFs) (last line = residual)
//              - NB_ITERATIONS: effective number of sifting iterations for each mode
//   Description
//
// cemdc_fix computes bivariate EMD, first algorithm [1] with NB_ITERATONS sifting
// iterations for each IMF
//
//   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
//
//   Bibliography
//
//
// [1] G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly.,
// "Bivariate Empirical Mode Decomposition",
// Signal Processing Letters (submitted)
//
//
//   Examples
//      X = rand(1,512)+%i*rand(1,512);
//      T=linspace(0,20,512);
//      
//     IMF = cemdc_fix(T,X);
//     [IMF,NB_IT] = cemdc_fix([],X);
//     IMF = cemdc_fix(T,X,20);
//     [IMF,NB_IT] = cemdc_fix([],X,[],4);
//
// See also
//  cemd_visu
//  emd_visu
//  emd
//  cemdc
//  cemdc2
//  cemdc2_fix
//
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007 gabriel.rilling@ens-lyon.fr
//
// code based on a student project by T. Boustane and G. Quellec, 11.03.2004
// supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
// email : pchainai@isima.fr).

