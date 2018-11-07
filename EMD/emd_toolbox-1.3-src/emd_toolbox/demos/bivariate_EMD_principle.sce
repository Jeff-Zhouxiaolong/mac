//bivariate_EMD_principle.m
//shows principle of the bivariate EMD extension
//reproduces Fig. 1 in "Bivariate Empirical Mode Decomposition", G. Rilling,
//P. Flandrin, P. Goncalves and J. M. Lilly, IEEE Signal Processing Letters
//
//  H. Nahrstaedt - Aug 2010
//G. Rilling 3/2007 email:  gabriel.rilling@ens-lyon.fr
//function bivariate_EMD_principle
N = 8192; // number of points
tmax = 40;
Ndirs1 = 2; // used to display the envelope curves
Ndirs2 = 8; // used to render a nice 3D tube envelope
stacksize('max')
if ~isdef("fmodany") then

function [y,iflaw]=fmodany(iflaw,t0);
// Signal with arbitrary frequency modulation.
// Calling Sequence
//	[Y,IFLAW]=fmodany(IFLAW)
//	[Y,IFLAW]=fmodany(IFLAW,T0)
//  Parameters
//	IFLAW : vector of the instantaneous frequency law samples.
//	T0    : time reference		(default: 1).
//	Y     : output signal
//    Description
//      fmodany generates a frequency modulated
//	signal whose instantaneous frequency law is approximately given by
//	the vector IFLAW (the integral is approximated by CUMSUM).
//	The phase of this modulation is such that y(t0)=1.
//    Examples
//       [y1,ifl1]=fmlin(100); [y2,ifl2]=fmsin(100);
//       iflaw=[ifl1;ifl2]; sig=fmodany(iflaw);
//       subplot(211); plot(real(sig))
//       subplot(212); plot(iflaw);
//   See also
//      fmhyp
//      fmsin
//      fmpar
//      fmconst
//      fmlin
//      fmpower
//   Authors
//      H. Nahrstaedt - Aug 2010
//	F. Auger, August 1995.
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
[ifrow,ifcol]=size(iflaw);
if (ifcol~=1),
error('IFLAW must have one column');
elseif (max(abs(iflaw))>0.5),
error('The elements of IFLAW should not be higher than 0.5');
end;
if (nargin==1),
t0=1;
elseif (t0==0)|(t0>ifrow),
error('T0 should be between 1 and length(IFLAW)');
end;

y=exp(%i*2.0*%pi*cumsum(iflaw));
y=y*conj(y(t0));
endfunction

end;
// display parameters
A = [0,tmax,-10,10,-10,10];
L = 2; // thickness of the envelope curves
gray = [.5,.5,.5];
target_faces_number = 1e3; // ~number of faces of the tube patch object. reduce/increase this number to get faster/nicer rendering
propNames = {
    'AmbientStrength'
    'DiffuseStrength'
    'SpecularStrength'
    'SpecularExponent'
    'SpecularColorReflectance'
    'FaceColor'
    'FaceAlpha'
    'EdgeColor'
  };
//propVals = {0,0,3,7,1,gray,.3,'none'};
Axprops = {
  'XTickLabel'
  'YTickLabel'
  'ZTickLabel'
  'CameraPosition'
  'CameraViewAngle'
  };
CameraPosition = [-69.1664 -166.3868 17.2944];
CameraViewAngle = 7.4163;
Axpropvals = {[],[],[],CameraPosition,CameraViewAngle};

function subplot_label(label)
label_pos = [.07,.07,.09];
tmax = 40;A = [0,tmax,-10,10,-10,10];
    xstring(A(1)+(A(2)-A(1))*label_pos(1),A(3)+(A(4)-A(3))*(1-label_pos(2)),A(5)+(A(6)-A(5))*(1-label_pos(3)),['('+label+')']);
endfunction
// abscissa vector
absc = linspace(0,tmax,N);
dt = mean(diff(absc));
// extended abscissa to avoid edge effects on the final plot
absc_calc = [-5:dt:-dt,absc,tmax+dt:dt:tmax+5];
tt = find(absc_calc>=0 & absc_calc <=tmax);

// fast rotating component
function y=fm1(t,tmax)
 y = 1+.4*sin(2*%pi/tmax*t);
endfunction
function y=am1(t,tmax)
y = 1+.5*cos(%pi/(2*tmax)*t);
endfunction
function y=dirstretch1(t,tmax)
y =  3*exp(%i*t*%pi/(1.5*tmax));
endfunction
x1 = emd_dirstretch(am1(absc_calc,tmax).*fmodany(dt*fm1(absc_calc,tmax).').',dirstretch1(absc_calc,tmax));

// slow rotating component
function y=fm2(t)
y =  0.0382+t-t;
endfunction
function y=am2(t)
y =  .2*(1+t);
endfunction

x2 = am2(absc_calc).*fmodany(dt*fm2(absc_calc).').';

// composite signal
x = x1+x2;

// compute the envelope curves
// for a large number of directions to get a smooth tube envelope

[env2] = cenvelope(x,Ndirs2);
// env2 = flipud([emax;emin]);
[env1,moy] = cenvelope(x,Ndirs1);
// env1 = flipud([emax;emin]);


function y = fliplr(x)
//  Copyright Aldo I Maalouf

if ndims(x)~=2,
disp('X must be a 2-D matrix!')
end
y = x(:,$:-1:1);
endfunction

// make MATLAB render an extended tube envelope in order to get
// good normals at the boundaries of the final patch object
[m,n] = size(env2);
env = [env2;env2(1:3,:)];
step = max(round(length(env)/target_faces_number),1);
inds = [fliplr(tt(1)-step:-step:1),tt(1):step:n];
if inds($)<n
    inds = [inds,n];
end
inds2 = find(inds>=tt(1) & inds<= tt($));
PYZ = env(:,inds);
PX = repmat(absc_calc(inds),m+3,1);
tmp = figure;
surf(PX,real(PYZ),imag(PYZ));
//VN = get(h,'VertexNormals');
//VN = VN(2:$-1,inds2,:);
PYZ = PYZ(2:$-1,inds2);
PX = PX(2:$-1,inds2);
//delete(h);
//close(tmp);


f=scf()//'defaultAxesNextplot','add','renderer','opengl','name','Principle of the Bivariate EMD');
f.figure_name='Principle of the Bivariate EMD';
// NOTE: openGL rendering is required for transparency (at least on GNU/Linux)

// function plot3c(varargin)
// param3d(varargin(1),real(varargin(2)),imag(varargin(2)),varargin(3:$));
// endfunction
subplot(221);
// plot the signal
plot3c(absc,x(tt))
xgrid
//axis(A)
//set(gca,Axprops,Axpropvals);
//subplot_label('a')

 subplot(222);
// plot the signal
plot3c(absc,x(tt))
// and the envelope curves
for k = 1:2*Ndirs1
    plot3c(absc,env1(k,tt))//,'k')//,'LineWidth',L)
e=gce() //the handle on the 3D polyline
e.foreground=color('red');
end
//set(gca,Axprops,Axpropvals);
xgrid
// plot the tube envelope
surf(PX,real(PYZ),imag(PYZ))//,gray,'VertexNormals',VN);
//set(h,propNames,propVals);
//axis(A)
//l1 = light('Position',[25,15,20]);
//l2 = light('Position',[-50,10,-18]);
//subplot_label('b')

subplot(223);
// plot the mean of the envelope curves
plot3c(absc,x(tt)-moy(tt))
// and the zero axis
plot3c(absc,zeros(1,length(absc)))//,'k')
e=gce() //the handle on the 3D polyline
e.foreground=color('red');
//axis(A)
//set(gca,Axprops,Axpropvals);
xlabel('Time')
xgrid
//subplot_label('c')

subplot(224);
plot3c(absc,moy(tt))
//axis(A)
//set(gca,Axprops,Axpropvals);
xlabel('Time')
xgrid
//subplot_label('d')

//hlink = linkprop(ax,{'CameraPosition','CameraUpVector'});
