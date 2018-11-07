function [IMF,NB_ITERATIONS]=emdc_fix(T,X,NB_ITERATONS,MAX_IMFS);
//  computes Empirical Mode Decomposition
//  Calling Sequence
//    [IMF,NB_ITERATIONS]=emdc_fix(T,X);
//    [IMF,NB_ITERATIONS]=emdc_fix([],X);
//    [IMF,NB_ITERATIONS]=emdc_fix(T,X,NB_ITERATONS);
//    [IMF,NB_ITERATIONS]=emdc_fix(T,X,NB_ITERATONS,MAX_IMFS);
// Parameters
// inputs:
//       - T: sampling times. If T=[], the signal is assumed uniformly sampled.1xN time instants
//       - X: analyzed signal, 1xN signal data
//       - NB_ITERATIONS: number of sifting iterations to be performed to  extract each IMF. If NB_ITERATIONS is empty or unspecified, 10 iterations   are performed by default.  Note: The effective number of sifting iterations might be less   than NB_ITERATIONS for the last modes if the sifting process has   to be stopped  because of a lack of extrema.
//       - MAX_IMFS: maximum number of IMFs to be extracted. If MAX_IMFS is  zero, empty or unspecified, the default behavior is to extract as  many IMFs as possible.
// outputs:
//              - IMF: intrinsic mode functions (IMFs) (last line = residual)
//             - NB_ITERATIONS: effective number of sifting iterations for each mode
//   Description
//
//   emdc_fix computes EMD according to [1] with NB_ITERATONS sifting iterations for each IMF
//
//   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
//   &
//   |#zeros-#extrema|<=1
//
//   Bibliography
//
// [1] N. E. Huang et al., "The empirical mode decomposition and the
// Hilbert spectrum for non-linear and non stationary time series analysis",
// Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
//
//   Examples
//      X = rand(1,512);
//      T=linspace(0,20,512);
//
//    IMF = emdc_fix(T,X);
//    [IMF,NB_IT] = emdc_fix([],X);
//    IMF = emdc_fix(T,X,100);
//   [IMF,NB_IT] = emdc_fix([],X,[],4);
//
// See also
//  emd_visu
//  emd
//  emdc
//  hhspectrum
//
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007 gabriel.rilling@ens-lyon.fr
//
// code based on a student project by T. Boustane and G. Quellec, 11.03.2004
// supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
// email : pchainai@isima.fr).
