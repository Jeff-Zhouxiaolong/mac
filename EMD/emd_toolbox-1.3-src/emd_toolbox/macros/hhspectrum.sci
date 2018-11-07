function [A,f,tt] = hhspectrum(t,x,l,aff)
//  compute Hilbert-Huang spectrum
// Calling Sequence
// [A,f,tt] = hhspectrum(x)
// [A,f,tt] = hhspectrum([],x)
// [A,f,tt] = hhspectrum(t,x)
// [A,f,tt] = hhspectrum(x,l)
// [A,f,tt] = hhspectrum([],x,l)
// [A,f,tt] = hhspectrum(t,x,l)
// [A,f,tt] = hhspectrum(x,l,aff)
// [A,f,tt] = hhspectrum([],x,aff)
// [A,f,tt] = hhspectrum(t,x,l,aff)
// Parameters
// inputs:
//   - x   : matrix with one signal per row
//   - t   : time instants
//   - l   : estimation parameter for emd_instfreq (integer >=1 (1:default))
//   - aff : if 1, displays the computation evolution
// outputs:
//   - A   : instantaneous amplitudes
//   - f   : instantaneous frequencies
//   - tt  : truncated time instants
// Description
//  need the Time-Frequency Toolbox (stftb)
//Examples
//
//      s = rand(1,512,'normal');
//      imf = emd(s);
//      [A,f,tt] = hhspectrum(imf(1:$-1,:));
//     [im,tt]=toimage(A,f);
//     disp_hhs(im);
//
//     s = rand(10,512,'normal');
//    [A,f,tt] = hhspectrum([],s,2,1);
//     [im,tt]=toimage(A,f);
//     disp_hhs(im);
//
// See also
//  emd
//  toimage
//  disp_hhs
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr

[nargout,nargin]=argn(0);

// if (nargin == 0),
//  error ( 'The number of parameters must be at least 1.' );
// end;
//
// if nargin < 2
//
//   t=1:size(x,2);
//
// end
//
// if nargin < 3
//
//   l=1;
//
// end
//
// if nargin < 4
//
//   aff = 0;
//
// end

select nargin
 case 1
   x=t;
   t=1:size(x,2);
   l=1;
   aff = 0;
 case 2
   if sum(length(x))==1
     aff = 0;
     l=x;
     x=t;
     t=1:size(x,2);
   else
      l=1;
      aff = 0;
   end;
 case 3
   if sum(length(x))==1
      aff=l;
      l=x;
      x=t;
       t=1:size(x,2);
   else
       aff = 0;
    end;
end
if isempty(t)
     t=1:size(x,2);
end;

if min(size(x)) == 1
	if size(x,2) == 1
		x = x';
		if nargin < 2
			t = 1:size(x,2);
		end
	end
	Nmodes = 1;
else
	Nmodes = size(x,1);
end

lt=length(t);

tt=t((l+1):(lt-l));

for i=1:Nmodes

  an(i,:)=hilbert(x(i,:)')';
  f(i,:)=emd_instfreq(an(i,:)',tt,l)';
  A=abs(an(:,l+1:$-l));

  if aff
	disprog(i,Nmodes,max(Nmodes,100))
  end

end
endfunction
