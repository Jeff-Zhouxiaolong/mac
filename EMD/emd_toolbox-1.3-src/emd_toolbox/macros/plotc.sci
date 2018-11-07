function obj = plotc(varargin)
//  plots a complex signal in 2D projection with variable projection angle.
//  Calling Sequence
//  plotc
// Description
//The angle can be modified at any time using the slider:
// slider at bottom: angle=0 the projection is the real part of the signal
// slider at 1/4: angle=pi/2 the projection is the imaginary part of the signal
// slider at top: angle=2pi same as 0
//
//Additionally, the user can choose to lock the axes or not when the angle
//is changed through the axes context menu. Note that the latter is disabled
//when tools from the figure toolbar (zoom,...) are selected.
//
//
//  multiple uses of plotc in the same axes are possible: it
//      produces only one slider that controls all the plots.
// Examples
//     s = rand(1,512,'normal')+%i*rand(1,512,'normal');
//     plotc(1:length(s),s);
//     s2 = rand(1,512,'normal')+%i*rand(1,512,'normal');
//     plotc(1:length(s2),s2,'r');
// See also
//  plot3c
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr

 [nargout,nargin]=argn(0);
slider_prop = 0.05;
slider_dec = 0.01;
lock_axes = 0;

// input processing
if nargin >=2 & sum(length(varargin(1)))==1 & typeof(varargin(1))=="handle"
  ax = varargin(1);
 f = ax.parent;
  varargin = varargin(2:$);
else
 // f=scf();
  ax = gca();
f = ax.parent;
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
    if (size(y,1)==1 | size(y,2)==1)
      x = 1:length(y);
      y = y(:);
    else
      x = 1:size(y,1);
    end
  case 2
    y = varargin(2);
    x = varargin(1);
end

//args = varargin(numeric_args+1:$);


// context menu and slider
// cmenu = get(ax,'UIContextMenu');
// if isempty(cmenu)
//   cmenu = uicontextmenu;
//   set(ax,'UIContextMenu',cmenu);
// end
// if ~emd_hastag(ax,'complexplot')
if isempty(emd_findtag('complexplotslider'))
//   lock_axes_menu_item = uimenu(cmenu,'Label','Lock axes','Callback',@togglelock);
     P = ax.axes_bounds;//get(ax,'Position');
    H=0.9;
    L=0.1;
//   slider = uicontrol('Style','Slider','Min',0,'Max',1,'Value',0,'SliderStep',[.01,.05],'Units','normalized','Position',[P(1)+P(3)+P(3)*slider_dec,P(2),P(3)*slider_prop,P(4)],'Callback',@slider_callback);
     slider = uicontrol(f, "units","normalized", "style", "slider","string","Update","tag","tagSlider","position",[P(1)+P(3)-P(3)*slider_dec,P(2),P(3)*slider_prop,P(4)],"Min",1,"Max",101,"SliderStep",[1 1],"Callback",...
"slider_plotc_callback()");
   emd_addtag(slider,'complexplotslider')
 else
   slider = emd_findtag('complexplotslider');
 end

// setup the complex plot
phi = get_axes_phase();

if (ax==[])
plot(x,real(exp(-%i*phi)*y),varargin(numeric_args+1:$));
else
plot(ax,x,(exp(-%i*phi)*y),varargin(numeric_args+1:$));
end
obj=gce();
emd_addtag(ax,'complexplot') // needs to be after the plot command in case the latter redefines the Tag property of the axes (default behavior)
emd_addtag(obj,'complex');
obj.children.user_data=y(:);
//arrayfun(@(x)setappdata(obj(x),'complex_data',y(:,x)),1:length(obj));
//set(obj,'UIContextMenu',cmenu);
obj

endfunction
//set(fig,'toolbar','figure');

   function show_phase(phi)
      drawlater();
//     Xl = get(ax,'Xlim');
//     Yl = get(ax,'Ylim');
       ax_ch = emd_findtag('complex');
       ax_ch=ax_ch.parent;
       minn=%inf;maxx=-%inf;
       for child = 1:length(ax_ch.children)
//       z = getappdata(child,'complex_data');
         z=ax_ch.children(child).children.user_data;
         ax_ch.children(child).children.data(:,2) = real(exp(-%i*phi)*z);
         minn=min([minn;real(exp(-%i*phi)*z)]);
         maxx=max([maxx;real(exp(-%i*phi)*z)]);
//       set(child,'YData',z);
       end
//     set(ax,'Xlim',Xl);
//     if ~lock_axes
//       set(ax,'YlimMode','auto');
//     else
//       set(ax,'Ylim',Yl);
//     end
     ax_ch.data_bounds(:,2)=[minn;maxx];
    drawnow();
   endfunction
//
   function slider_plotc_callback()
     show_phase(get_axes_phase());
   endfunction
//
   function phi = get_axes_phase()
      //F_Slider = findobj('tag','tagSlider');
     F_Slider=  emd_findtag('complexplotslider');
      v = get(F_Slider,'Value');
//     phi = 2*pi*v;
       phi = 2*%pi*(v-1)/100;
   endfunction
//
//   function togglelock(varargin)
//     lock_axes = ~lock_axes;
//     if lock_axes
//       set(lock_axes_menu_item,'Checked','on');
//     else
//       set(lock_axes_menu_item,'Checked','off');
//     end
//   endfunction
