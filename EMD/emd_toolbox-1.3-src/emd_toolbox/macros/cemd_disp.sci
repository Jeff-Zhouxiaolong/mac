function cemd_disp(t,x,env,rmode)
// displays complex envelope curves and the corresponding tube envelope
//  Calling Sequence
//      cemd_disp(X,ENV)
//      cemd_disp([],X,ENV)
//      cemd_disp(T,X,ENV)
//      cemd_disp(X,ENV,MODE)
//      cemd_disp(T,X,ENV,MODE)
//      cemd_disp([],X,ENV,MODE)
// Parameters
// inputs:
//       - T: time instants
//       - X: analyzed signal (complex)
//       - ENV: matrix returned by cenvelope.m Each line is an envelope curve
//       - MODE: 'render' -> 3D rendering of the tube enclosing the signal    'wire'   -> wireframe display (much faster if there is no hardware acceleration)
//
// Description
// TODO: the program uses zbuffer rendering by default. You may switch to openGL
// rendering by commenting/uncommenting a line at the beginning of the function.
// openGL rendering is faster but zbuffer is generally nicer.
// Examples
//     s = rand(1,512,'normal')+%i*rand(1,512,'normal');
//     [env] = cenvelope(s);
//     cemd_disp(s,env);
// See also
//  cenvelope
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007 gabriel.rilling@ens-lyon.fr


 [nargout,nargin]=argn(0);
renderer = 'zbuffer'; // nicer
// renderer = 'openGL'; // faster

DEF_mode = 'render';

if nargin == 2
    env = x;
    x = t;
    t = 1:length(x);
end

if nargin == 3
    if type(env)==10
        rmode = env;
        env = x;
        x = t;
        t = 1:length(x);
    end
end
if isempty(t)
  t = 1:length(x);
end

t=t(:)';
x=x(:)';

if ~exists('rmode')
    rmode = DEF_mode;
end

//if (rmode=='render')
    target_faces_number = 5000;
//else
//    target_faces_number = 1000;
//end

//col = [.5,.5,.5];
h=scf();
h.figure_name='complex envelope';
plot3c(t,x)
xgrid

[m,n] = size(env);

step = max(round(length(env)/target_faces_number),1);
env = [env;env(1:3,:)];
inds = 1:step:n;
if inds($)<n
    inds = [inds,n];
end
PYZ = env(:,inds);
PX = mtlb_repmat(t(inds),m+3,1);

surf(PX,real(PYZ),imag(PYZ),'facecol','red','edgecol','blu");//,col);
//VN = get(h,'VertexNormals');
//VN = VN(2:$-1,:,:);
//PYZ = PYZ(2:$-1,:);
//PX = PX(2:$-1,:);
//delete(h);
//surf(PX,real(PYZ),imag(PYZ),col,'VertexNormals',VN);
//hidden off

//select convstr(part(rmode,(1:min(4,length(rmode)))))
//    case 'rend'
        //set(h,'EdgeColor','none')
       // set(h,'FaceColor',col)
        for k = 1:size(env(1:$-3,:),1)
            plot3c(t,env(k,:))
        end
//    case 'wire'
        //set(h,'EdgeColor',col)
       // set(h,'FaceColor','none')
//    else
//        warning('cemd_disp:input unknown option '+rmode)
//end
//set(h,'FaceLighting','phong')
//set(gcf,'renderer',renderer)

//l(1) = light('Position',[min(PX(:)),min(real(PYZ(:))),min(imag(PYZ(:)))]);
//l(2) = light('Position',[max(PX(:)),max(real(PYZ(:))),max(imag(PYZ(:)))]);

// if nargout > 0
//     varargout = {h,l};
// else
//     varargout = {};
// end
endfunction

