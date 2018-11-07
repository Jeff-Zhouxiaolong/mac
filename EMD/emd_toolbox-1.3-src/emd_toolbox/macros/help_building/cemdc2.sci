function [IMF,NB_ITERATIONS]=cemdc2(T,X,STOP_PARAMETERS,MAX_IMFS,NDIRS)
//  bivariate Empirical Mode Decomposition, second algorithm
//  Calling Sequence
// [IMF,NB_ITERATIONS]=cemdc2(T,X);
// [IMF,NB_ITERATIONS]=cemdc2([],X);
// [IMF,NB_ITERATIONS]=cemdc2(T,X,STOP_PARAMETERS);
// [IMF,NB_ITERATIONS]=cemdc2(T,X,STOP_PARAMETERS,MAX_IMFS);
// [IMF,NB_ITERATIONS]=cemdc2(T,X,STOP_PARAMETERS,MAX_IMFS,NDIRS);
// Parameters
// inputs:
//       - T: sampling times. If T=[], the signal is assumed uniformly sampled.1xN time instants
//       - X: analyzed signal,   1xN signal data 
//       - STOP_PARAMETERS: parameters for the stopping criterion:  if scalar the value is used to specify THRESHOLD only.  otherwise the vector should be: [THRESHOLD,TOLERANCE].   if STOP_PARAMETERS is unspecified or empty, default values are used: [0.05,0.05]
//       - MAX_IMFS: maximum number of IMFs to be extracted. If MAX_IMFS is    zero, empty or unspecified, the default behavior is to extract as   many IMFs as possible.
//       - NDIRS: number of directions used to compute the local mean.   If unspecified, the default value is 4.  TODO: the actual number of directions (according to [1]) is 2*NDIRS
// outputs: 
//              - IMF: intrinsic mode functions (IMFs) (last line = residual)
//              - NB_ITERATIONS: effective number of sifting iterations for each mode
//   Description
// cemdc2 computes bivariate EMD, second algorithm [1] with stopping criterion for
// sifting similar to the one proposed in [2]:
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
// [2] G. Rilling, P. Flandrin and P. Gonçalves
// "On Empirical Mode Decomposition and its algorithms",
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003
//
//   Examples
//      X = rand(1,512)+%i*rand(1,512);
//      T=linspace(0,20,512);
//      
//   IMF = cemdc2(T,X);
//   [IMF,NB_IT] = cemdc2([],X);
//   IMF = cemdc2(T,X,0.1);
//   IMF = cemdc2(T,X,[0.1,0.1]);
//   [IMF,NB_IT] = cemdc2([],X,[],4);
//
//
// See also
//   emd_visu
//   cemd_visu
//   emd
//   cemdc2_fix
//   cemdc
//   cemdc_fix 
//
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007 gabriel.rilling@ens-lyon.fr
//
// code based on a student project by T. Boustane and G. Quellec, 11.03.2004
// supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
// email : pchainai@isima.fr).
