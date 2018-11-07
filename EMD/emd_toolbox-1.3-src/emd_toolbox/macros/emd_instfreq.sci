function [fnormhat,t]=emd_instfreq(x,t,L,ptrace);
// Instantaneous frequency estimation.
// Calling Sequence
//	[FNORMHAT,T]=emd_instfreq(X)
//	[FNORMHAT,T]=emd_instfreq(X,T)
//	[FNORMHAT,T]=emd_instfreq(X,T,L)
//	[FNORMHAT,T]=emd_instfreq(X,T,L,TRACE)
//  Parameters
//	X : Analytic signal to be analyzed.
//	T : Time instants	        (default : 2:length(X)-1).
//	L : If L=1, computes the (normalized) instantaneous frequency  of the signal X defined as moc_angle(X(T+1)*conj(X(T-1)) ;   if L>1, computes a Maximum Likelihood estimation of the   instantaneous frequency of the deterministic part of the signal    blurried in a white gaussian noise.	    L must be an integer       	(default : 1).
//	TRACE : if nonzero, the progression of the algorithm is shown   (default : 0).
//	FNORMHAT : Output (normalized) instantaneous frequency.
//	T : Time instants.
//  Description
//      emd_instfreq computes the instantaneous frequency of the analytic signal X at time instant(s) T, using the
//	trapezoidal integration rule.
//	The result FNORMHAT lies between 0.0 and 0.5.
//   Authors
//      H. Nahrstaedt - Aug 2010
//	F. Auger, March 1994, July 1995.
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
error('At least one parameter required');
end;
[xrow,xcol] = size(x);
if (xcol~=1),
error('X must have only one column');
end

if (nargin == 1),
t=2:xrow-1; L=1; ptrace=0.0;
elseif (nargin == 2),
L = 1; ptrace=0.0;
elseif (nargin == 3),
ptrace=0.0;
end;

if L<1,
error('L must be >=1');
end
[trow,tcol] = size(t);
if (trow~=1),
error('T must have only one row');
end;

if (L==1),
if or(t==1)|or(t==xrow),
error('T can not be equal to 1 neither to the last element of X');
else
fnormhat=0.5*(moc_angle(-x(t+1).*conj(x(t-1)))+%pi)/(2*%pi);
end;
else
H=kaytth(L);
if or(t<=L)|or(t+L>xrow),
error('The relation L<T<=length(X)-L must be satisfied');
else
for icol=1:tcol,
//if ptrace, disprog(icol,tcol,10); end;
ti = t(icol); tau = 0:L;
R = x(ti+tau).*conj(x(ti-tau));
M4 = R(2:L+1).*conj(R(1:L));

d=2e-6;
tetapred = H * (moc_unwrap(moc_angle(-M4))+%pi);
while tetapred<0.0 , tetapred=tetapred+(2*%pi); end;
while tetapred>2*%pi, tetapred=tetapred-(2*%pi); end;
iter = 1;
while (d > 1e-6)&(iter<50),
M4bis=M4 .* exp(-%i*2.0*tetapred);
teta = H * (moc_unwrap(moc_angle(M4bis))+2.0*tetapred);
while teta<0.0 , teta=(2*%pi)+teta; end;
while teta>2*%pi, teta=teta-(2*%pi); end;
d=abs(teta-tetapred);
tetapred=teta; iter=iter+1;
end;
fnormhat(icol,1)=teta/(2*%pi);
end;
end;
end;

endfunction

function H=kaytth(flength);
//	 Kay-Tretter filter computation.
// Calling Sequence
//	H=kaytth(length)
//    See also
//        emd_instfreq
//   Authors
//      H. Nahrstaedt - Aug 2010
//	F. Auger, March 1994.
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

pp1=flength*(flength+1);
den=2.0*flength*(flength+1)*(2.0*flength+1.0)/3.0;
i=1:flength; H=pp1-i.*(i-1);

H=H ./ den;

endfunction

function disprog(i,N,steps);
// Display progression of a loop.
// Calling Sequence
//	disprog(i,N,steps)
// Parameters
//	I     : loop variable
//	N     : final value of i
//	STEPS : number of displayed steps.
// Description
//       disprog displays the progression of a loop.
// Examples
//      N=16; for i=1:N, disprog(i,N,5); end;
// Authors
//      H. Nahrstaedt - Aug 2010
//	F. Auger, August, December 1995.
//       from an idea of R. Settineri.
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

if (i==1),
tic();
end;

if (i==N),
printf('100 %% complete in %g seconds.\n', toc())
elseif (floor(i*steps/N)~=floor((i-1)*steps/N)),
printf('%g ', floor(i*steps/N)*ceil(100/steps));
end;

endfunction
