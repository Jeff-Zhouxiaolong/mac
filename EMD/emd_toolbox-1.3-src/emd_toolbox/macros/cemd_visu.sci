function  cemd_visu(t,x,imf,plottype,i)
//  visualization of bivariate/complex EMD
// Calling Sequence
//    cemd_visu(x,imf)
//    cemd_visu([],x,imf)
//    cemd_visu(t,x,imf)
//    cemd_visu(x,imf,plottype)
//    cemd_visu([],x,imf,plottype)
//    cemd_visu(t,x,imf,plottype)
//    cemd_visu(x,imf,plottype,i)
//    cemd_visu([],x,imf,plottype,i)
//    cemd_visu(t,x,imf,plottype,i)
// Parameters
// inputs :
//            - x : complex analyzed signal
//            - t : time instants
//            - imf : set of complex IMFs
//            - plottype (optional) :
// 			'imf' - plots imf
// 			'c2f' - plots c2f
// 			'f2c' - plots f2c
//            - i (optional) : figure number for display

// Description
//
//  emd_reconstruct can be used to calculate fine to coarse and coarse to fine reconstruction
//
// The slider on the right hand side of the figure controls the direction on
// which the complex signal is projected:
// phi=0: the projection corresponds to the real part of the signal,
// phi=pi/2: the projection corresponds to the imaginary part of the signal.
//
// Examples
//     s = rand(1,512,'normal')+%i*rand(1,512,'normal');
//     imf = cemdc([],s);
//     cemd_visu(1:length(s),s,imf);
//
// See also
//   emd_reconstruct
//  emd
//  cemdc
//  cemdc_fix
//   cemdc2
//   cemdc2_fix
// Authors
// H. Nahrstaedt
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr

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
        plottype='imf';
  end
 case 5
    fignum = i;
end
if isempty(t)
    t=1:length(x);
end;
clear i;
plottype=convstr(plottype);

if plottype=="all"
 plottype="imf";
end
if isempty(t)
  t = 1:length(x);
end

// if nargin == 2
//     imf = x;
//     x = t;
//     t = 1:length(x);
//      h=scf();
// end
//
// if nargin == 3
//     if sum(length(imf)==1
//         numfig = imf;
//         imf = x;
//         x = t;
//         t = 1:length(x);
//         h=scf(numfig)
//      else
//        h=scf();
//     end
// else
//   if(nargin==4)
//     h=scf(numfig)
//   else
//     h=scf()
//   end
//
// end


if length(x)~= length(t)
  error('X and T must have the same length');
end
if size(imf,2) ~= length(x)
  error('the number of columns in IMF must equal the length of X')
end


// [k,n] = size(imf);
//
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
phi = 0;

if plottype=="c2f"
    imf=c2f;
elseif plottype=="f2c"
    imf=f2c;
end;

if exists('fignum')
  h=scf(fignum);clf()
else
  h=scf();clf()
end
//set(gcf,'name','cemd_visu: EMD','Toolbar','figure')
if plottype=="imf"
h.figure_name='cemd_visu: EMD';
elseif plottype=="f2c"
h.figure_name=('f2c : fine to coarse reconstruction ')
elseif plottype=="c2f"
h.figure_name=('c2f : coarse to fine reconstruction')
end
h.tag = "cemd_visu";


mx = mean(x);
Mx = max(abs(x-mx));
ind=find(~isnan(imf($,:)));
mr=mean(imf($,ind));
Mr = max(abs(imf($,ind)-mr));
M = max(max(abs(imf(1:k-1,:))));
A_x = [t(1),t(n),real(mx)-Mx,real(mx)+Mx];
A_imf = [t(1),t(n),-M,M];
A_res = [t(1),t(n),real(mr)-Mr,real(mr)+Mr];
ax_imf=list();
subplot(k+1,1,1);
ax_x=gca();
//set(gca,'NextPlot','replacechildren','YTick',[],'XTickLabel',{},'Box','on')
//axis(A_x);
ax_x.data_bounds=A_x;
ylabel(['signal'])
title('Empirical Mode Decomposition, \phi/2\pi=0');

for j=1:k-1
  subplot(k+1,1,j+1);
  ax_imf(j) = gca();
  //set(gca,'NextPlot','replacechildren','YTick',[],'XTickLabel',{},'Box','on')
   ylabel([plottype+string(j)]);
  //axis(A_imf);
  ax_imf(j).data_bounds=A_imf;
end
subplot(k+1,1,k+1);
 ax_res = gca();
//set(gca,'NextPlot','replacechildren','Box','on')
//axis(A_res);
ax_res.data_bounds=A_res;
if plottype=="imf"
  ylabel('res.')
else
   ylabel("signal");
end

plot_emd(ax_x,ax_imf,ax_res,x,t,imf,phi);

P(1,:) = ax_x.axes_bounds;
P(2,:) = ax_res.axes_bounds;
L = (P(1,3)-P(1,1));
L= 0.01;
H = P(1,2)+P(1,4)-P(2,2);
H=0.9;
L=0.1;

global handles;
handles.dummy = 0;
handles.ax_x=ax_x;
handles.ax_imf=ax_imf;
handles.ax_res=ax_res;
handles.x=x;
handles.t=t;
handles.imf=imf;

//slider_handle = uicontrol('Style','Slider','Min',0,'Max',1,'Value',0,'SliderStep',[.01,.05],'Units','normalized','Position',[1-2*L/3,P(2,2),L/3,H],'Callback',@slider_callback);
handles.slider = uicontrol(h, "units","normalized", "style", "slider","string","Update","tag","tagSlider","position",[1-2*L/3,P(2,1),L/3,H],"Min",1,"Max",101,"SliderStep",[1 1],"Callback",...
"slider_callback()");

 // varargout = list(f2c,c2f);
endfunction
   function plot_emd(ax_x,ax_imf,ax_res,x,t,imf,phi)
    drawlater();
    if (length(ax_x.children)>0),
       ax_x.children.children.data(:,2)=real(exp(-%i*phi)*x(:));
    else
       plot(ax_x,t,real(exp(-%i*phi)*x));
    end;
    //set(ax_x,'Ylim',[real(exp(-i*phi)*mx)-Mx,real(exp(-i*phi)*mx)+Mx]);
    //ax_x.data_bounds(:,2)=[real(exp(-%i*phi)*mx)-Mx;real(exp(-%i*phi)*mx)+Mx];
    [k,n] = size(imf);
    for j = 1:k-1
      if (length(ax_imf(j).children)>0),
	ax_imf(j).children.children.data(:,2)=real(exp(-%i*phi)*imf(j,:))';
      else
	plot(ax_imf(j),t,real(exp(-%i*phi)*imf(j,:)));
      end
    end
    if (length(ax_res.children)>0),
      ax_res.children.children.data(:,2)=real(exp(-%i*phi)*imf($,:))';
    else
      plot(ax_res,t,real(exp(-%i*phi)*imf($,:)),'r');
    end
    drawnow();
     //h = findobj('tag','cemd_visu');
     //h.figure_name=['Empirical Mode Decomposition, \phi/2\pi='+string(phi/(2*%pi))];
    //set(ax_res,'Ylim',[real(exp(-i*phi)*mr)-Mr,real(exp(-i*phi)*mr)+Mr]);
     //ax_res.data_bounds(:,2)=[real(exp(-%i*phi)*mr)-Mr;real(exp(-%i*phi)*mr)+Mr];
    //set(T,'String',['Empirical Mode Decomposition, \phi/2\pi=',num2str(phi/(2*pi))]);
    //drawnow;
  endfunction


   function slider_callback()
      global handles;
//     v = get(slider_handle,'Value');
      F_Slider = findobj('tag','tagSlider');
      v = get(F_Slider,'Value')
      //v=handles.slider.Value;
      //phi = 2*%pi*(handles.slider.Value-1)/100;
      phi = 2*%pi*(v-1)/100;
     plot_emd(handles.ax_x,handles.ax_imf,handles.ax_res,handles.x,handles.t,handles.imf,phi);
   endfunction
