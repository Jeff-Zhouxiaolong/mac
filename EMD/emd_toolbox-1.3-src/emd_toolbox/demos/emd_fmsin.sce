// P. Flandrin, Mar. 13, 2003 - modified Mar. 2, 2006
//  H. Nahrstaedt - Aug 2010
//
clc;  h=scf(); fig_id=get(gcf(),"figure_id"); clf; mode(1);lines(0);
// computes and displays EMD for the sum of 2 sinusoidal
// FM's + 1 Gaussian logon
//
// displays reassigned spectrograms of the sum signal and of
// the 3 first modes extracted
//
// produces Figure 1 in
//
// G. Rilling, P. Flandrin and P. Gon�alv�s
// "On Empirical Mode Decomposition and its algorithms"
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003
warning("This demo needs the Time Frequency Toolbox!");
mode(-1)

if ~isdef("fmsin") then
function [y,iflaw]=fmsin(N,fnormin,fnormax,period,t0,fnorm0,pm1);
//	Signal with sinusoidal frequency modulation.
// Calling Sequence
//	[Y,IFLAW]=fmsin(N)
//	[Y,IFLAW]=fmsin(N,FNORMIN)
//	[Y,IFLAW]=fmsin(N,FNORMIN,FNORMAX)
//	[Y,IFLAW]=fmsin(N,FNORMIN,FNORMAX,PERIOD)
//	[Y,IFLAW]=fmsin(N,FNORMIN,FNORMAX,PERIOD,T0)
//	[Y,IFLAW]=fmsin(N,FNORMIN,FNORMAX,PERIOD,T0,FNORM0)
//	[Y,IFLAW]=fmsin(N,FNORMIN,FNORMAX,PERIOD,T0,FNORM0,PM1)
// Parameters
//	N       : number of points.
//	FNORMIN : smallest normalized frequency          (default: 0.05)
//	FNORMAX : highest normalized frequency           (default: 0.45)
//	PERIOD  : period of the sinusoidal fm            (default: N   )
//	T0      : time reference for the phase           (default: N/2 )
//	FNORM0  : normalized frequency at time T0        (default: 0.25)
//	PM1     : frequency direction at T0 (-1 or +1)	 (default: +1  )
//	Y       : signal
//	IFLAW   : its instantaneous frequency law (optional).
//  Description
//    fmsin generates a frequency modulation with a sinusoidal frequency.
//	This sinusoidal modulation is designed such that the instantaneous
//	frequency at time T0 is equal to FNORM0, and the ambiguity
//	between increasing or decreasing frequency is solved by PM1.
//   Examples
//      z=fmsin(140,0.05,0.45,100,20,0.3,-1.0);plot(real(z));
//   See also
//      fmhyp
//      fmpower
//      fmodany
//      fmconst
//      fmlin
//      fmpar
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
error ( 'At least 1 parameter required' ) ;
elseif (nargin == 1),
fnormin = 0.05; fnormax = 0.45; period = N;
t0= round(N/2); fnorm0 = 0.25; pm1=+1;
elseif (nargin == 2),
fnormax = 0.45; period = N;
t0 = round(N/2); fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 3),
period = N; t0 = round(N/2);
fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 4),
t0 = round(N/2); fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 5),
fnorm0 = 0.5*(fnormin+fnormax); pm1=+1;
elseif (nargin == 6),
pm1= +1;
elseif (nargin==7),
if (abs(pm1)~=1),
error('pm1 must be equal to -1 or +1');
end;
end;

if (N <= 0),
error ('The signal length N must be strictly positive' );
elseif (abs(fnormin) > 0.5)|(abs(fnormax) > 0.5)|(abs(fnorm0) > 0.5),
error ('fnormin, fnormax and fnorm0 must be between -0.5 and 0.5') ;
elseif (fnormin > fnormax)
error ('fnormin must be lower than fnormax');
elseif (fnormin > fnorm0)|(fnorm0 > fnormax),
error ('fnorm0 must be between fnormin and fnormax') ;
else
fnormid=0.5*(fnormax+fnormin);
delta  =0.5*(fnormax-fnormin);
phi=-pm1*acos((fnorm0-fnormid)/delta);
time=(1:N)-t0;
phase=2*%pi*fnormid*time+delta*period*(sin(2*%pi*time/period+phi)-sin(phi));
y = exp(%i*phase).';
if (nargout==2)
iflaw=fnormid+delta*cos(2*%pi*time'/period+phi);
end
end
endfunction

end;
if ~isdef("amgauss")
function y = amgauss(N,t0,T);
// Generate gaussian amplitude modulation.
// Calling Sequence
//	Y=amgauss(N)
//	Y=amgauss(N,T0)
//	Y=amgauss(N,T0,T)
// Parameters
//	N  : number of points.
//	T0 : time center		(default : N/2).
//	T  : time spreading		(default : 2*sqrt(N)).
//	Y  : signal.
// Description
//  amgauss  generates a gaussian amplitude modulation
//	centered on a time T0, and with a spread proportional to T.
//	This modulation is scaled such that Y(T0)=1
//	and Y(T0+T/2) and Y(T0-T/2) are approximately equal to 0.5 .
// Examples
//       z=amgauss(160); plot(z);
//       z=amgauss(160,90,40); plot(z);
//       z=amgauss(160,180,50); plot(z);
//  See also
//      amexpo1s
//      amexpo2s
//      amrect
//      amtriang
// Authors
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
t0=N/2; T=2*sqrt(N);
elseif (nargin ==2),
T=2*sqrt(N);
end;

if (N<=0),
error('N must be greater or equal to 1.');
else
tmt0=(1:N)'-t0;
y = exp(-(tmt0/T).^2 * %pi);
end;
endfunction
end;
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

N = 2000;// # of data samples
T = 1:4:N;
t = 1:N;

p = N/2;// period of the 2 sinusoidal FM's

// sinusoidal FM 1
fmin1 = 1/64;// min frequency
fmax1 = 1.5*1/8;// max frequency
x1 = fmsin(N,fmin1,fmax1,p,N/2,fmax1);

// sinusoidal FM 1
fmin2 = 1/32;// min frequency
fmax2 = 1.5*1/4;// max frequency
x2 = fmsin(N,fmin2,fmax2,p,N/2,fmax2);

// logon
f0 = 1.5*1/16;// center frequency
x3 = amgauss(N,N/2,N/8).*fmconst(N,f0);

a1 = 1;
a2 = 1;
a3 = 1;

x = real(a1*x1+a2*x2+a3*x3);
x = x/max(abs(x));

[imf,ort,nbits] = emd(x);

emd_visu(t,x,imf,'all',fig_id);


// time-frequency distributions
Nf = 256;// # of frequency bins
Nh = 127;// short-time window length
w = tftb_window(Nh,'Kaiser');

[s,rs] = tfrrsp(x,T,Nf,w,1);
[s,rs1] = tfrrsp(imf(1,:)',T,Nf,w,1);
[s,rs2] = tfrrsp(imf(2,:)',T,Nf,w,1);
[s,rs3] = tfrrsp(imf(3,:)',T,Nf,w,1);

scf(fig_id+3);clf();
f = scf();f.color_map = jetcolormap(64);
subplot(221)
grayplot(1:500,1:128,rs(1:128,:)')
//set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('signal')



subplot(222)
grayplot(1:500,1:128,rs1(1:128,:)')

//set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #1')



subplot(223)
grayplot(1:500,1:128,rs2(1:128,:)')
//set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #2')


subplot(224)
grayplot(1:500,1:128,rs3(1:128,:)')
//set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('time')
ylabel('frequency')
title('mode #3')
