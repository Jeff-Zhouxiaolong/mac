function [im,tt,ff] = toimage(A,f,param1,param2,param3)
//  transforms a spectrum made of 1D functions in an 2D image
// Calling Sequence
//[im,tt,ff] = toimage(A,f)
//[im,tt,ff] = toimage(A,f,t)
//[im,tt,ff] = toimage(A,f,splx)
//[im,tt,ff] = toimage(A,f,splx,sply)
//[im,tt,ff] = toimage(A,f,t,splx)
//[im,tt,ff] = toimage(A,f,t,splx,sply)
// Parameters
// inputs :
//            - A    : amplitudes of modes (1 mode per row of A)
//            - f    : instantaneous frequencies
//            - t    : time instants
//            - splx : number of columns of the output im (time resolution).If different from length(t), works only for uniform  sampling.
//            - sply : number of rows of the output im (frequency resolution).
// outputs :
//            - im   : 2D image of the spectrum
//            - tt   : time instants in the image
//            - ff   : centers of the frequency bins
// Description
// toimage transforms a spectrum made of 1D functions (e.g., output of "hhspectrum") in an 2D image
// Examples
//     s = rand(1,512,'normal');
//     imf = emd(s);
//     [A,f,tt] = hhspectrum(imf(1:$-1,:));
//     [im,tt]=toimage(A,f);
//     disp_hhs(im);
// See also
//  emd
//  hhspectrum
//  disp_hhs
// Authors
// H. Nahrstaedt - Aug 2010 - 2013
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr


[nargout,nargin]=argn(0);

if (nargin < 2),
 error ( 'The number of parameters must be at least 2.' );
end;

DEFSPL = 400;



select nargin
  case 2
    t = 1:size(A,2);
    sply = DEFSPL;
    splx = length(t);
  case 3
    if (sum(length(param1))==1)
      t = 1:size(A,2);
      splx = length(t);
      sply = param1;
    else
      t = param1;
      splx = length(t);
      sply = DEFSPL;
    end
  case 4
    if (sum(length(param1))==1)
      t = 1:size(A,2);
      sply = param1;
      splx = param2;
    else
      t = param1;
      sply = param2;
      splx = length(t);
    end
  case 5
    t = param1;
    splx = param2;
    sply = param3;
end
if (size(A,1)==1 | size(A,2)==1)
  A = A(:)';
  f = f(:)';
end


if type(A)==5 | ~isreal(A) | length(size(A)) > 2
  error('A argument must be a real matrix')
end
if type(f)==5 | ~isreal(f) | length(size(f)) > 2
  error('f argument must be a real matrix')
end
if or(size(f)~=size(A))
  error('A and f matrices must have the same size')
end
if type(t)==5 | ~isreal(t) | ~(size(t,1)==1 | size(t,2)==1) | length(t)~=size(A,2)
  error('t argument must be a vector and its length must be the number of columns in A and f inputs')
end
if ~(sum(length(splx))==1) | ~isreal(splx) | splx ~= floor(splx) | splx <= 0
  error('splx argument must be a positive integer')
end
if ~(sum(length(sply))==1) | ~isreal(sply) | sply ~= floor(sply) | sply <= 0
  error('splx argument must be a positive integer')
end

if or(diff(diff(t))) & splx ~= length(t)
  warning('toimage:nonuniformtimeinsants','When splx differs from length(t), the function only works for equally spaced time instants. You may consider reformating your data (using e.g. interpolation) before using toimage.')
end

f = min(f,0.5);
f = max(f,0);

indf = round(2*f*(sply-1)+1);
indt = repmat(round(linspace(1,length(t),splx)),size(A,1),1);
im = moc_accumarray([indf(:),indt(:)],A(:),[sply,splx]);

indt = indt(1,:);
tt = t(indt);
ff = (0:sply-1)*0.5/sply+1/(4*sply);

endfunction
