function emd_visu(t,x,imf,plottype,i);
//  visualization of EMD and partial reconstructions (fine to coarse and coarse to fine)
// Calling Sequence
//    emd_visu(x,imf)
//    emd_visu([],x,imf)
//    emd_visu(t,x,imf)
//    emd_visu(x,imf,plottype)
//    emd_visu([],x,imf,plottype)
//    emd_visu(t,x,imf,plottype)
//    emd_visu(x,imf,plottype,i)
//    emd_visu([],x,imf,plottype,i)
//   emd_visu(t,x,imf,plottype,i)
// Parameters
// inputs :
//           - x : analyzed signal, if x is complex cemd_visu is called
//            - t : time instants
//            - imf : output of emd.m
//            - plottype (optional) :
// 			'all' (default) - plots imf, c2f and f2c
// 			'imf' - plots imf
// 			'c2f' - plots c2f
// 			'f2c' - plots f2c
//            - i (optional) : figure number for display
// Description
//  emd_reconstruct can be used to calculate fine to coarse and coarse to fine reconstruction
// Examples
//     s = rand(1,512,'normal');
//     imf = emd(s);
//     emd_visu(s,imf);
// See also
//   emd_reconstruct
//  cemd_visu
//  emd
//  emdc
//  emdc_fix
//  cemdc
//  cemdc_fix
//  cemdc2
//  cemdc2_fix
// Authors
// H. Nahrstaedt - Aug  2010
// P. Flandrin, Mar. 13, 2003
// G. Rilling, last modification 3.2006 gabriel.rilling@ens-lyon.fr

[nargout,nargin]=argn(0);

select nargin
 case 1
  error("At least 2 parameters needed!");
 case 2
   imf=x;
   x=t;
   t=1:length(x);
   plottype='all';
 case 3
   if type(imf)==10
     plottype=imf;
     imf=x;
     x=t;
      t=1:length(x);
  else
   plottype='all';
  end;
 case 4
   if type(imf)==10
    i=plottype;
    plottype=imf;
    imf=x;
    x=t;
    t=1:length(x);
     fignum = i;
    elseif type(plottype)~=10
      i=plottype;
       fignum = i;
       plottype='all';
  end
 case 5
    fignum = i;
end
if isempty(t)
    t=1:length(x);
end;
clear i;
plottype=convstr(plottype);
// if sum(size(t)>1)>1
//   imf = t;
//   t = 1:length(x);
//   if(nargin==3)
//     fignum = i;
//     plottype='all';
//   end
// else
//   if(nargin==4)
//     fignum = i;
//   end
// end

if ~isreal(x)|~isreal(imf)
  if exists('fignum')
    cemd_visu(t,x,imf,plottype,fignum);
  else
    cemd_visu(t,x,imf,plottype);
  end
  return;
end

if length(x)~= length(t)
  error('X and T must have the same length');
end
if size(imf,2) ~= length(x)
  error('the number of columns in IMF must equal the length of X')
end

if plottype=="all" | plottype=="imf" then
if exists('fignum')
  scf(fignum);clf()
else
  scf();clf()
end

mx = min(x);
Mx = max(x);

s = size(imf);
k = s(1);

M = max(max(abs(imf(1:k-1,:))));

subplot(k+1,1,1)
plot(t,x)
a=gca();
a.data_bounds=[t(1) t(s(2)) mx Mx];
//axis([t(1) t(s(2)) mx Mx])
//set(gca,'YTick',[])
//set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,imf(j,:))
  a=gca();a.data_bounds=[t(1) t(s(2)) -M M];
  //axis([t(1) t(s(2)) -M M])
  //set(gca,'YTick',[])
  //set(gca,'XTick',[])
  ylabel(['imf'+string(j)])
end
subplot(k+1,1,1)
title('Empirical Mode Decomposition')

subplot(k+1,1,k+1)
plot(t,imf(k,:),'r')
//axis('tight')
//set(gca,'YTick',[])
//set(gca,'XTick',[])
ylabel('res.')
end
// f2c = [];
// f2c(1,:) = imf(1,:);
//
// c2f = [];
// c2f(1,:) = imf(k,:);
//
// for j = 2:k
//   f2c(j,:) = f2c(j-1,:) + imf(j,:);
//   c2f(j,:) = c2f(j-1,:) + imf(k+1-j,:);
// end
[f2c, c2f] = emd_reconstruct(imf);

if plottype=="all" | plottype=="f2c" then
if exists('fignum')
  scf(fignum+1)
else
  scf()
end

mx = min(x);
Mx = max(x);

s = size(f2c);
k = s(1);

M = max(max(abs(f2c(1:k-1,:))));

subplot(k+1,1,1)
plot(t,x)
 a=gca();a.data_bounds=[t(1) t(s(2)) mx Mx];
//axis([t(1) t(s(2)) mx Mx])
//set(gca,'YTick',[])
//set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,f2c(j,:))
  //axis([t(1) t(s(2)) -M M])
  a=gca();a.data_bounds=[t(1) t(s(2)) -M M];
  //set(gca,'YTick',[])
  //set(gca,'XTick',[])
  ylabel(['f2c'+string(j)])
end
subplot(k+1,1,1)
title('f2c : fine to coarse reconstruction ')

subplot(k+1,1,k+1)
plot(t,f2c(k,:),'r')
mr = min(f2c(k,:));
Mr = max(f2c(k,:));
  a=gca();a.data_bounds=[t(1) t(s(2)) mr Mr];
//axis([t(1) t(s(2)) mr Mr])
//set(gca,'YTick',[])
//set(gca,'XTick',[])
ylabel('signal')
end

if plottype=="all" | plottype=="c2f" then
if exists('fignum')
  scf(fignum+2)
else
  scf();
end

mx = min(x);
Mx = max(x);

s = size(c2f);
k = s(1);

M = max(max(abs(c2f(1:k-1,:)-mean(x))));

subplot(k+1,1,1)
plot(t,x)
a=gca();a.data_bounds=[t(1) t(s(2)) mx Mx];
// axis([t(1) t(s(2)) mx Mx])
// set(gca,'YTick',[])
// set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,c2f(j,:))
  a=gca();a.data_bounds=[t(1) t(s(2)) -M M];
//   axis([t(1) t(s(2)) -M M])
//   set(gca,'YTick',[])
//   set(gca,'XTick',[])
  ylabel(['c2f'+string(j)])
end
subplot(k+1,1,1)
title('c2f : coarse to fine reconstruction')

subplot(k+1,1,k+1)
plot(t,c2f(k,:),'r')
mr = min(c2f(k,:));
Mr = max(c2f(k,:));
  a=gca();a.data_bounds=[t(1) t(s(2)) mr Mr];
// axis([t(1) t(s(2)) mr Mr])
// set(gca,'YTick',[])
// set(gca,'XTick',[])
ylabel('signal')
end
//varargout = list(f2c,c2f);
endfunction
