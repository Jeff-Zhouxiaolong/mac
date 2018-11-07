function disp_hhs(im,param1,param2,param3)
//  display Hilbert-Huang spectrum
// Calling Sequence
//      disp_hhs(im)
//      disp_hhs(im,t)
//      disp_hhs(im,threshold)
//      disp_hhs(im,t,threshold)
//      disp_hhs(im,threshold,fs)
//      disp_hhs(im,[],fs)
//      disp_hhs(im,t,[],fs)
//      disp_hhs(im,t,threshold,fs)
//      disp_hhs(im,t,threshold)
// Parameters
//          - im: image matrix (e.g., output of "toimage")
//          - t (optional): time instants (e.g., output of "toimage") 
//          - threshold: is the visualization threshold, in %   default: threshold=5;
//          - fs: sampling frequency
// Desciption
// displays in a new figure the spectrum contained in matrix "im" (amplitudes in dB).
// Examples 
//     s = rand(1,512,'normal');
//     imf = emd(s);
//     [A,f,tt] = hhspectrum(imf(1:$-1,:));
//     [im,tt]=toimage(A,f);
//     disp_hhs(im);
// See also
//  emd
//  hhspectrum
//  toimage
// Authors
// H. Nahrstaedt
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr


[nargout,nargin]=argn(0);

if (nargin == 0),
 error ( 'The number of parameters must be at least 1.' );
end;
fs = 0;
threshold = 5;

t = 1:size(im,2);
select nargin
  case 1
    
  case 2
    if (sum(length(param1))==1)
      threshold = param1;
    else
      t = param1;
    end
  case 3
    if ~(sum(length(param1))==1)
      t = param1;
      threshold = param2;
    else
      threshold = param1;
      fs = param2;
    end
  case 4
    t = param1;
    threshold = param2;
    fs = param3;
end



if threshold < 0
  threshold = -threshold;
elseif threshold == 0 | threshold==[]
  error('threshold must be nonzero')
end
threshold=min(100,threshold);
maxx = max(im);
 minn=max(min(im),maxx*threshold/100.0);


indmin=find(im<minn);
  im(indmin)=minn*ones(1,length(indmin));
 
  indmax=find(im>maxx);
  im(indmax)=maxx*ones(1,length(indmax));
// warning off

 // im(find(im==0))=-1e10;
 // im(find(im~=-1e10)) = 10*log10(im(find(im~=-1e10))/M);

// warning on

f=scf(); 
f.color_map = jetcolormap(32);
if fs == 0
  //colorbar(inf,0);
  grayplot(t,linspace(0,0.5,size(im,1)),log10(im)');//,[inf,0]);
  ylabel('normalized frequency')
else
  //colorbar(inf,0);
  grayplot(t,linspace(0,0.5*fs,size(im,1)),log10(im)');//,[inf,0]);
  ylabel('frequency')
end
//set(gca,'YDir','normal')
xlabel('time')
title('Hilbert-Huang spectrum')
endfunction