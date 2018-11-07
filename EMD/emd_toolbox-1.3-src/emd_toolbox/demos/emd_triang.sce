// EMD_TRIANG.M
//
// P. Flandrin, Mar. 13, 2003 - modified Mar. 2, 2006
//  H. Nahrstaedt - Aug 2010
//
clc;  h=scf(); fig_id=get(gcf(),"figure_id"); clf; mode(1);lines(0);
// computes and displays EMD for the sum of 2 triangular
// waveforms + 1 tone
//
// produces Figure 2 in
//
// G. Rilling, P. Flandrin and P. Gon�alv�s
// "On Empirical Mode Decomposition and its algorithms"
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003
mode(-1)

N = 1024;// # of data samples
t = 1:N;
demopath = get_absolute_file_path("emd_triang.sce");
exec(demopath+"triangular_signal.sci",-1);
if ~isdef("fmconst") then
function [y,iflaw] = fmconst(N,fnorm,t0);
// Signal with constant frequency modulation.
// Calling Sequence
//	[Y,IFLAW] = fmconst(N,FNORM,T0)
//     Parameters
//	N     : number of points.
//	FNORM : normalised frequency.       (default: 0.25)
//	T0    : time center.                (default: round(N/2))
//	Y     : signal.
//	IFLAW : instantaneous frequency law (optional).
//     Description
//      fmconst generates a frequency modulation  with a constant frequency fnorm.
//	The phase of this modulation is such that y(t0)=1.
//     Examples
//        z=amgauss(128,50,30).*fmconst(128,0.05,50); plot(real(z));
//
//     See also
//      fmlin
//      fmsin
//      fmodany
//      fmhyp
//      fmpar
//      fmpower
//   Authors
//      H. Nahrstaedt - Aug 2010
//	F. Auger, July 1995.
//	Copyright (c) 1996 by CNRS (France).

//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
[nargout,nargin]=argn(0);
if (nargin == 0),
error ( 'The number of parameters must be at least 1.' );
elseif (nargin == 1),
t0=round(N/2); fnorm=0.25 ;
elseif (nargin ==2),
t0=round(N/2); // this has been changed from N/2 to round N/2 for matlab 5
end;

if (N<=0),
error('N must be greater or equal to 1.');
elseif (abs(fnorm)>0.5),
error('The normalised frequency must be between -0.5 and 0.5');
else
tmt0=(1:N)'-t0;
y = exp(%i*2.0*%pi*fnorm*tmt0);
y=y/y(t0);
if (nargout==2), iflaw=fnorm*ones(N,1); end;
end;

endfunction
end;


// triangular waveform 1
p1 = fix(N/6);// period
x1 = triangular_signal(N,p1);

// tone
f0 = 0.03;// frequency
x2 = real(fmconst(N,f0));

// triangular waveform 2
p3 = 5;// period
x3 = triangular_signal(N,p3);

x = x1 + 0.4*x2' + .4*x3;

[imf,ort,nbits] = emd(x);

emd_visu(t,x,imf,'all',fig_id);
