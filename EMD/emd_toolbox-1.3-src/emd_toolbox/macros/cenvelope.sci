function [env,envmoy] = cenvelope(t,x,Nphases,INTERP)
//  computes envelope curves for bivariate EMD
// Calling Sequence
//      [env, mean] = cenvelope(x);
//      [env, mean] = cenvelope([],x);
//      [env, mean] = cenvelope(x,ndirs);
//      [env, mean] = cenvelope([],x,ndirs);
//      [env, mean] = cenvelope(t,x,ndirs);
//      [env, mean] = cenvelope(x,ndirs,INTERP);
//      [env, mean] = cenvelope([],x,ndirs,INTERP);
//      [env, mean] = cenvelope(t,x,ndirs,INTERP);
// Parameters
// inputs :
//          - x: analyzed signal
//          - t (optional): sampling times, default 1:length(x)
//          - ndirs: number of directions used to compute the mean (default: 4)   TODO: the actual number of directions according to the paper is 2*ndirs
//          - interp (optional): interpolation scheme: 'linear', 'cubic' or 'spline' (default)
// outputs :
//           - env: each stands for an envelope curve sustaining the tube envelope of    the complex signal
//           - mean = mean of the envelope curves (corresponding to the first algorithm in the paper)
//
// Description
//  cenvelope computes envelope curves for bivariate EMD [1]
// Bibliography
// [1] G. Rilling, P. Flandrin, P. Gon�alves and J. M. Lilly.,
// "Bivariate Empirical Mode Decomposition",
// Signal Processing Letters (submitted)
//
// See also
//  cemd_disp
//  emd
//  cemdc
//  cemdc_fix
//  cemdc2
//  cemdc2_fix
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr

[nargout,nargin]=argn(0);

NBSYM = 2;
DEF_INTERP = 'spline';

select nargin
 case 1
  x=t;
  t=1:length(x);
  Nphases = 4;
  INTERP = DEF_INTERP;
 case 2
   if sum(length(x))==1
      Nphases=x;
      x=t;
      t=1:length(x);
     INTERP = DEF_INTERP;
  else
    Nphases = 4;
    INTERP = DEF_INTERP;
  end
 case 3
    if sum(length(x))==1
      INTERP=Nphases;
      Nphases=x;
      x=t;
      t=1:length(x);
  else
    INTERP = DEF_INTERP;
  end
end

if isempty(t)
    t=1:length(x);
end;



// if nargin < 2
//     x = t;
//     t = 1:length(x);
// end
//
// if nargin >= 2
//     if (sum(length(x))==1)
//         if nargin == 3
//           INTERP = Nphases;
//         end
//         Nphases = x;
//         x = t;
//         t = 1:length(x);
//     end
// end

if ~exists('INTERP')
    INTERP = DEF_INTERP;
end
if ~exists('Nphases')
    Nphases = 4;
end
if isempty(t)
     t=1:length(x);
end;

for k = 1:Nphases
    phi = (k-1)*%pi/Nphases;
    y = real(exp(-%i*phi)*x);
    //[im,iM] = extr(y);
    [im,iM] = emd_local_peaks(y);
    [tmin,tmax,zmin,zmax] = boundary_conditions_emd(im(:)',iM(:)',t(:)',y(:)',x(:)',NBSYM);
    envmin(k,:) = interp1(real(tmin),real(zmin),t,INTERP);
    envmax(k,:) = interp1(real(tmax),real(zmax),t,INTERP);
end

env=[envmin;envmax];

if nargout == 2
  envmoy = mean(envmax + envmin,1)/2;
end

endfunction
