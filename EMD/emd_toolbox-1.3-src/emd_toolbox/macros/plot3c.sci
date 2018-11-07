function plot3c(varargin)
//  plots a complex signal in 3D
//  Calling Sequence
//  plot3c
// Description
// use: same as param3d
// See also
//  plotc
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification 12.2006 gabriel.rilling@ens-lyon.fr

 [nargout,nargin]=argn(0);
// if first arument is a handle to some axes object
if nargin >=2 & sum(length(varargin(1)))==1 & typeof(varargin(1))=="handle"
  ax = varargin(1);
  varargin = varargin(2:$);
else
  ax = [];//gca();
end

// char_inds = find(cellfun(@ischar,varargin));
// if ~isempty(char_inds)
//   numeric_args = char_inds(1)-1;
// else
//   numeric_args = length(varargin);
// end
numeric_args=0;
for k=1:length(varargin)
  if type(varargin(k))==1
    numeric_args=numeric_args+1;
  end
end

select numeric_args
  case 1
    y = varargin(1);
    x = 1:length(y);
  case 2
    y = varargin(2);
    x = varargin(1);
end
y=y(:)';
x=x(:)';
if ax==[]
  param3d(x,real(y),imag(y),varargin(numeric_args+1:$));
else
 param3d(ax,x,real(y),imag(y),varargin(numeric_args+1:$));
end
endfunction