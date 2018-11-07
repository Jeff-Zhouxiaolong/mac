function [imf,ort,nbits] = emd(x,varargin)
//  computes Empirical Mode Decomposition
// Calling Sequence
// [IMF,ORT,NB_ITERATIONS] = emd(X)
// [IMF,ORT,NB_ITERATIONS] = emd(X,...,'Option_name',Option_value,...)
// [IMF,ORT,NB_ITERATIONS] = emd(X,OPTS)
//   Parameters
//  stopping criterion options:
//
// STOP: vector of stopping parameters [THRESHOLD,THRESHOLD2,TOLERANCE] if the input vector's length is less than 3, only the first parameters are set, the remaining ones taking default values. default: [0.05,0.5,0.05]
//
// FIX (int): disable the default stopping criterion and do exactly <FIX>  number of sifting iterations for each mode
//
// FIX_H (int): disable the default stopping criterion and do <FIX_H> sifting  iterations with |#zeros-#extrema|<=1 to stop [4]
//
//  bivariate/complex EMD options:
//
// COMPLEX_VERSION: selects the algorithm used for complex EMD ([3])
// COMPLEX_VERSION = 1: "algorithm 1"
// COMPLEX_VERSION = 2: "algorithm 2" (default)
//
// NDIRS: number of directions in which envelopes are computed (default 4) rem: the actual number of directions (according to [3]) is 2*NDIRS
//
//  other options:
//
// T: sampling times (line vector) (default: 1:length(x))
//
// MAXITERATIONS: maximum number of sifting iterations for the computation of each mode (default: 2000)
//
// MAXMODES: maximum number of imfs extracted (default: Inf)
//
// DISPLAY: if equals to 1 shows sifting steps with pause if equals to 2 shows sifting steps without pause (movie style) rem: display is disabled when the input is complex
//
// INTERP: interpolation scheme: 'linear', 'cubic', 'pchip' or 'spline' (default) see interp1 documentation for details
//
// MASK: masking signal used to improve the decomposition according to [5]
//
// Description
//
// IMF = EMD(X) where X is a real vector computes the Empirical Mode Decomposition [1] of X, resulting in a matrix IMF containing 1 IMF per row, the/ last one being the residue. The default stopping criterion is the one proposed in [2]:
//
//   at each point, mean_amplitude < THRESHOLD2*envelope_amplitude
//   &
//   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
//   &
//   |#zeros-#extrema|<=1
//
// where mean_amplitude = abs(envelope_max+envelope_min)/2
// and envelope_amplitude = abs(envelope_max-envelope_min)/2
//
// IMF = EMD(X) where X is a complex vector computes Bivariate Empirical Mode
// Decomposition [3] of X, resulting in a matrix IMF containing 1 IMF per row, the
// last one being the residue. The default stopping criterion is similar to the
// one proposed in [2]:
//
//   at each point, mean_amplitude < THRESHOLD2*envelope_amplitude
//   &
//   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
//
// where mean_amplitude and envelope_amplitude have definitions similar to the
// real case
//
// IMF = EMD(X,...,'Option_name',Option_value,...) sets options Option_name to
// the specified Option_value (see Options)
//
// IMF = EMD(X,OPTS) is equivalent to the above syntax provided OPTS is a struct
// object with field names corresponding to option names and field values being the
// associated values
//
// [IMF,ORT,NB_ITERATIONS] = EMD(...) returns an index of orthogonality
//
// and the number of iterations to extract each mode in NB_ITERATIONS
//   Bibliography
//
//
// [1] N. E. Huang et al., "The empirical mode decomposition and the
// Hilbert spectrum for non-linear and non stationary time series analysis",
// Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
//
// [2] G. Rilling, P. Flandrin and P. Gonzalves
// "On Empirical Mode Decomposition and its algorithms",
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003
//
// [3] G. Rilling, P. Flandrin, P. Gonzalves and J. M. Lilly.,
// "Bivariate Empirical Mode Decomposition",
// Signal Processing Letters (submitted)
//
// [4] N. E. Huang et al., "A confidence limit for the Empirical Mode
// Decomposition and Hilbert spectral analysis",
// Proc. Royal Soc. London A, Vol. 459, pp. 2317-2345, 2003
//
// [5] R. Deering and J. F. Kaiser, "The use of a masking signal to improve
// empirical mode decomposition", ICASSP 2005
//   Examples
//X = rand(1,512);
//
//IMF = emd(X);
//
//IMF = emd(X,'STOP',[0.1,0.5,0.05],'MAXITERATIONS',100);
//
//T=linspace(0,20,1e3);
//X = 2*exp(%i*T)+exp(3*%i*T)+.5*T;
//IMF = emd(X,'T',T);
// cemd_visu(X,IMF);
//
//OPTIONS.DISPLAY = 1;
//OPTIONS.FIX = 10;
//OPTIONS.MAXMODES = 3;
//[IMF,ORT,NBITS] = emd(real(X),OPTIONS);
//
//
//
// See also
//  emd_visu
//  hhspectrum
//  emdc
//  emdc_fix
//  cemdc
//  cemdc_fix
//  cemdc2
//  cemdc2_fix
// Authors
// Holger Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007 gabriel.rilling@ens-lyon.fr

    [nargout,nargin]=argn(0);
if (nargin==0)
  error("At least one parameter needed!");
end

[x,t,sd,sd2,tol,MODE_COMPLEX,ndirs,display_sifting,sdt,sd2t,r,imf,k,nbit,NbIt,MAXITERATIONS,FIXE,FIXE_H,MAXMODES,INTERP,mask] = init(x,varargin);

if display_sifting
  fig_h = figure();
end


//main loop : requires at least 3 extrema to proceed
while ~stop_EMD(r,MODE_COMPLEX,ndirs) & (k < MAXMODES+1 | MAXMODES == 0) & ~or(mask)

  // current mode
  m = r;

  // mode at previous iteration
  mp = m;

  //computation of mean and stopping criterion
  if FIXE
    [stop_sift,moyenne] = stop_sifting_fixe(t,m,INTERP,MODE_COMPLEX,ndirs);
  elseif FIXE_H
    stop_count = 0;
    [stop_sift,moyenne] = stop_sifting_fixe_h(t,m,INTERP,stop_count,FIXE_H,MODE_COMPLEX,ndirs);
  else
    [stop_sift,moyenne] = stop_sifting(m,t,sd,sd2,tol,INTERP,MODE_COMPLEX,ndirs);
  end

  // in case the current mode is so small that machine precision can cause
  // spurious extrema to appear
  if (max(abs(m))) < (1e-10)*(max(abs(x)))
    if ~stop_sift
      warning('emd: '+'forced stop of EMD : too small amplitude')
    else
      disp('forced stop of EMD : too small amplitude')
    end
    break
  end


  // sifting loop
  while ~stop_sift & nbit<MAXITERATIONS

    if(~MODE_COMPLEX & nbit>MAXITERATIONS/5 & pmodulo(nbit,floor(MAXITERATIONS/10))==0 & ~FIXE & nbit > 100)
      disp(['mode '+string(k)+', iteration '+string(nbit)])
      if exists('s')
        disp(['stop parameter mean value : '+string(s)])
      end
      //[im,iM] = extr(m);
      [im,iM]=emd_local_peaks(m);
      disp(string(sum(m(im) > 0))+' minima > 0; '+string(sum(m(iM) < 0))+' maxima < 0.')
    end

    //sifting
    m = m - moyenne;

    //computation of mean and stopping criterion
    if FIXE
      [stop_sift,moyenne] = stop_sifting_fixe(t,m,INTERP,MODE_COMPLEX,ndirs);
    elseif FIXE_H
      [stop_sift,moyenne,stop_count] = stop_sifting_fixe_h(t,m,INTERP,stop_count,FIXE_H,MODE_COMPLEX,ndirs);
    else
      [stop_sift,moyenne,s] = stop_sifting(m,t,sd,sd2,tol,INTERP,MODE_COMPLEX,ndirs);
    end

    // display
    if display_sifting & ~MODE_COMPLEX
      NBSYM = 2;
      //[indmin,indmax] = extr(mp);
      [indmin,indmax] = emd_local_peaks(mp);
      [tmin,tmax,mmin,mmax] = boundary_conditions_emd(indmin,indmax,t,mp,mp,NBSYM);
      envminp = interp1(tmin,mmin,t,INTERP);
      envmaxp = interp1(tmax,mmax,t,INTERP);
      envmoyp = (envminp+envmaxp)/2;
      if FIXE | FIXE_H
        display_emd_fixe(t,m,mp,r,envminp,envmaxp,envmoyp,nbit,k,display_sifting,fig_h)
      else
        sxp=2*(abs(envmoyp))./(abs(envmaxp-envminp));
        sp = mean(sxp,'m');
        display_emd(t,m,mp,r,envminp,envmaxp,envmoyp,s,sp,sxp,sdt,sd2t,nbit,k,display_sifting,stop_sift,fig_h)
      end
    end

    mp = m;
    nbit=nbit+1;
    NbIt=NbIt+1;

    if(nbit==(MAXITERATIONS-1) & ~FIXE & nbit > 100)
      if exists('s')
        warning('emd '+['forced stop of sifting : too many iterations... mode '+string(k)+'.',' stop parameter mean value : '+string(s)])
      else
        warning('emd '+['forced stop of sifting : too many iterations... mode '+string(k)+'.'])
      end
    end

  end // sifting loop
  imf(k,:) = m;
  if display_sifting
    disp(['mode '+string(k)+' stored'])
  end
  nbits(k) = nbit;
  k = k+1;


  r = r - m;
  nbit=0;


end //main loop

if or(r) & ~or(mask)
  imf(k,:) = r;
end

ort = emd_io(x,imf);

if display_sifting
  close(fig_h);
end
endfunction

//---------------------------------------------------------------------------------------------------
// tests if there are enough (3) extrema to continue the decomposition
function stop = stop_EMD(r,MODE_COMPLEX,ndirs)
if MODE_COMPLEX
  for k = 1:ndirs
    phi = (k-1)*%pi/ndirs;
    //[indmin,indmax] = extr(real(exp(%i*phi)*r));
    [indmin,indmax] = emd_local_peaks(real(exp(%i*phi)*r));
    ner(k) = length(indmin) + length(indmax);
  end
  stop = or(ner < 3);
else
  //[indmin,indmax] = extr(r);
  [indmin,indmax] = emd_local_peaks(r);
  ner = length(indmin) + length(indmax);
  stop = ner < 3;
end
endfunction

//---------------------------------------------------------------------------------------------------
// computes the mean of the envelopes and the mode amplitude estimate
function [envmoy,nem,nzm,amp] = mean_and_amplitude(m,t,INTERP,MODE_COMPLEX,ndirs)
  [nargout,nargin]=argn(0);
NBSYM = 2;
if MODE_COMPLEX
  select MODE_COMPLEX
    case 1
      for k = 1:ndirs
        phi = (k-1)*%pi/ndirs;
        y = real(exp(-%i*phi)*m);
        //[indmin,indmax,indzer] = extr(y);
        [indmin,indmax] = emd_local_peaks(y);
        [indzer] = emd_zero_crossings(y);
        nem(k) = length(indmin)+length(indmax);
        nzm(k) = length(indzer);
        [tmin,tmax,zmin,zmax] = boundary_conditions_emd(indmin,indmax,t,y,m,NBSYM);
        envmin(k,:) = interp1(tmin,zmin,t,INTERP);
        envmax(k,:) = interp1(tmax,zmax,t,INTERP);
      end
      envmoy = mean((envmin+envmax)/2,1);
      if nargout > 3
        amp = mean(abs(envmax-envmin),1)/2;
      end
    case 2
      for k = 1:ndirs
        phi = (k-1)*%pi/ndirs;
        y = real(exp(-%i*phi)*m);
        //[indmin,indmax,indzer] = extr(y);
        [indmin,indmax] = emd_local_peaks(y);
        [indzer] = emd_zero_crossings(y);
        nem(k) = length(indmin)+length(indmax);
        nzm(k) = length(indzer);
        [tmin,tmax,zmin,zmax] = boundary_conditions_emd(indmin,indmax,t,y,y,NBSYM);
        envmin(k,:) = exp(%i*phi)*interp1(tmin,zmin,t,INTERP);
        envmax(k,:) = exp(%i*phi)*interp1(tmax,zmax,t,INTERP);
      end
      envmoy = mean((envmin+envmax),1);
      if nargout > 3
        amp = mean(abs(envmax-envmin),1)/2;
      end
  end
else
  //[indmin,indmax,indzer] = extr(m);
  [indmin,indmax] = emd_local_peaks(m);
  [indzer] = emd_zero_crossings(m);
  nem = length(indmin)+length(indmax);
  nzm = length(indzer);
  [tmin,tmax,mmin,mmax] = boundary_conditions_emd(indmin,indmax,t,m,m,NBSYM);
  envmin = interp1(tmin,mmin,t,INTERP);
  envmax = interp1(tmax,mmax,t,INTERP);
  envmoy = (envmin+envmax)/2;
  if nargout > 3
    amp = mean(abs(envmax-envmin),1)/2;
  end
end
endfunction

//-------------------------------------------------------------------------------
// default stopping criterion
function [stop,envmoy,s] = stop_sifting(m,t,sd,sd2,tol,INTERP,MODE_COMPLEX,ndirs)
try
  [envmoy,nem,nzm,amp] = mean_and_amplitude(m,t,INTERP,MODE_COMPLEX,ndirs);
  sx = abs(envmoy)./amp;
  s = mean(sx,'m');
  stop = ~(((sum(sx > sd)/length(sx)) > tol | or(sx > sd2)) & (and(nem > 2)));
  if ~MODE_COMPLEX
    stop = stop & ~(abs(nzm-nem)>1);
  end
  stop=sum(stop);
catch
  stop = 1;
  envmoy = zeros(1,length(m));
  s = %nan;
end
endfunction

//-------------------------------------------------------------------------------
// stopping criterion corresponding to option FIX
function [stop,moyenne]= stop_sifting_fixe(t,m,INTERP,MODE_COMPLEX,ndirs)
try
  moyenne = mean_and_amplitude(m,t,INTERP,MODE_COMPLEX,ndirs);
  stop = 0;
catch
  moyenne = zeros(1,length(m));
  stop = 1;
end
endfunction

//-------------------------------------------------------------------------------
// stopping criterion corresponding to option FIX_H
function [stop,moyenne,stop_count]= stop_sifting_fixe_h(t,m,INTERP,stop_count,FIXE_H,MODE_COMPLEX,ndirs)
try
  [moyenne,nem,nzm] = mean_and_amplitude(m,t,INTERP,MODE_COMPLEX,ndirs);
  if (all(abs(nzm-nem)>1))
    stop = 0;
    stop_count = 0;
  else
    stop_count = stop_count+1;
    stop = (stop_count == FIXE_H);
  end
catch
  moyenne = zeros(1,length(m));
  stop = 1;
end
endfunction

//-------------------------------------------------------------------------------
// displays the progression of the decomposition with the default stopping criterion
function display_emd(t,m,mp,r,envmin,envmax,envmoy,s,sb,sx,sdt,sd2t,nbit,k,display_sifting,stop_sift,fig_h)
scf(fig_h);clf(fig_h);
subplot(4,1,1)
plot(t,mp);//hold on
plot(t,envmax,'--k');plot(t,envmin,'--k');plot(t,envmoy,'r');
title(['IMF '+string(k)+';   iteration '+string(nbit)+' before sifting']);
//set(gca,'XTick',[])
//hold  off
subplot(4,1,2)
plot(t,sx)
//hold on
plot(t,sdt,'--r')
plot(t,sd2t,':k')
title('stop parameter')
//set(gca,'XTick',[])
//hold off
subplot(4,1,3)
plot(t,m)
title(['IMF '+string(k)+';   iteration '+string(nbit)+' after sifting']);
//set(gca,'XTick',[])
subplot(4,1,4);
plot(t,r-m)
title('residue');
disp(['stop parameter mean value : '+string(sb)+' before sifting and '+string(s)+' after'])
if stop_sift
  disp('last iteration for this mode')
end
if display_sifting == 2
  //pause(0.01)
  //pause
else
  //pause
end
endfunction

//---------------------------------------------------------------------------------------------------
// displays the progression of the decomposition with the FIX and FIX_H stopping criteria
function display_emd_fixe(t,m,mp,r,envmin,envmax,envmoy,nbit,k,display_sifting,fig_h)
scf(fig_h);clf(fig_h)
subplot(3,1,1)
plot(t,mp);//hold on;
plot(t,envmax,'--k');plot(t,envmin,'--k');plot(t,envmoy,'r');
title(['IMF '+string(k)+';   iteration ',+string(nbit)+' before sifting']);
//set(gca,'XTick',[])
//hold  off
subplot(3,1,2)
plot(t,m)
title(['IMF '+string(k)+';   iteration '+string(nbit)+' after sifting']);
//set(gca,'XTick',[])
subplot(3,1,3);
plot(t,r-m)
title('residue');
if display_sifting == 2
  //pause(0.01)
  //pause
else
  //pause
end
endfunction

//---------------------------------------------------------------------------------------
//defines new extrema points to extend the interpolations at the edges of the
//signal (mainly mirror symmetry)
// function [tmin,tmax,zmin,zmax] = boundary_conditions(indmin,indmax,t,x,z,nbsym)
//
// 	lx = length(x);
//
// 	if (length(indmin) + length(indmax) < 3)
// 		error('not enough extrema')
// 	end
//
//     // boundary conditions for interpolations :
//
//
// 	if indmax(1) < indmin(1)
// 		if x(1) > x(indmin(1))
// 			lmax = moc_fliplr(indmax(2:min(length(indmax),nbsym+1)));
// 			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
// 			lsym = indmax(1);
// 		else
// 			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
// 			lmin = [moc_fliplr(indmin(1:min(length(indmin),nbsym-1))),1];
// 			lsym = 1;
// 		end
// 	else
//
// 		if x(1) < x(indmax(1))
// 			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
// 			lmin = moc_fliplr(indmin(2:min(length(indmin),nbsym+1)));
// 			lsym = indmin(1);
// 		else
// 			lmax = [moc_fliplr(indmax(1:min(length(indmax),nbsym-1))),1];
// 			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
// 			lsym = 1;
// 		end
// 	end
//
// 	if indmax($) < indmin($)
// 		if x($) < x(indmax($))
// 			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
// 			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym,1):$-1));
// 			rsym = indmin($);
// 		else
// 			rmax = [lx,moc_fliplr(indmax(max(length(indmax)-nbsym+2,1):$))];
// 			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
// 			rsym = lx;
// 		end
// 	else
// 		if x($) > x(indmin($))
// 			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym,1):$-1));
// 			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
// 			rsym = indmax($);
// 		else
// 			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
// 			rmin = [lx,moc_fliplr(indmin(max(length(indmin)-nbsym+2,1):$))];
// 			rsym = lx;
// 		end
// 	end
// 	tlmin = 2*t(lsym)-t(lmin);
// 	tlmax = 2*t(lsym)-t(lmax);
// 	trmin = 2*t(rsym)-t(rmin);
// 	trmax = 2*t(rsym)-t(rmax);
//
// 	// in case symmetrized parts do not extend enough
// 	if tlmin(1) > t(1) | tlmax(1) > t(1)
// 		if lsym == indmax(1)
// 			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
// 		else
// 			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
// 		end
// 		if lsym == 1
// 			error('bug')
// 		end
// 		lsym = 1;
// 		tlmin = 2*t(lsym)-t(lmin);
// 		tlmax = 2*t(lsym)-t(lmax);
// 	end
//
// 	if trmin($) < t(lx) | trmax($) < t(lx)
// 		if rsym == indmax($)
// 			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
// 		else
// 			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
// 		end
// 		if rsym == lx
// 			error('bug')
// 		end
// 		rsym = lx;
// 		trmin = 2*t(rsym)-t(rmin);
// 		trmax = 2*t(rsym)-t(rmax);
// 	end
// 	zlmax =z(lmax);
// 	zlmin =z(lmin);
// 	zrmax =z(rmax);
// 	zrmin =z(rmin);
//
// 	tmin = [tlmin t(indmin) trmin];
// 	tmax = [tlmax t(indmax) trmax];
// 	zmin = [zlmin z(indmin) zrmin];
// 	zmax = [zlmax z(indmax) zrmax];
// endfunction
//
//
//
// function y = moc_fliplr(x)
// //  Copyright Aldo I Maalouf
//
// if ndims(x)~=2,
// disp('X must be a 2-D matrix!')
// end
// y = x(:,$:-1:1);
// endfunction

//---------------------------------------------------------------------------------------------------
//extracts the indices of extrema
// function [indmin, indmax, indzer] = extr(x,t)
// [nargout,nargin]=argn(0);
// if(nargin==1)
//   t=1:length(x);
// end
//
// m = length(x);
//
// if nargout > 2
//   x1=x(1:m-1);
//   x2=x(2:m);
//   indzer = find(x1.*x2<0);
//
//   if or(x == 0)
//     iz = find( x==0 );
//     indz = [];
//     if or(diff(iz)==1)
//       zer = x == 0;
//       dz = diff([0 zer 0]);
//       debz = find(dz == 1);
//       finz = find(dz == -1)-1;
//       indz = round((debz+finz)/2);
//     else
//       indz = iz;
//     end
//     indzer = mtlb_sort([indzer indz]);
//   end
// end
//
// d = diff(x);
//
// n = length(d);
// d1 = d(1:n-1);
// d2 = d(2:n);
// indmin = find(d1.*d2<0 & d1<0)+1;
// indmax = find(d1.*d2<0 & d1>0)+1;
//
//
// // when two or more successive points have the same value we consider only one extremum in the middle of the constant area
// // (only works if the signal is uniformly sampled)
//
// if or(d==0)
//
//   imax = [];
//   imin = [];
//
//   bad = (d==0);
//   dd = diff([0 bad 0]);
//   debs = find(dd == 1);
//   fins = find(dd == -1);
//   if debs(1) == 1
//     if length(debs) > 1
//       debs = debs(2:$);
//       fins = fins(2:$);
//     else
//       debs = [];
//       fins = [];
//     end
//   end
//   if length(debs) > 0
//     if fins($) == m
//       if length(debs) > 1
//         debs = debs(1:($-1));
//         fins = fins(1:($-1));
//
//       else
//         debs = [];
//         fins = [];
//       end
//     end
//   end
//   lc = length(debs);
//   if lc > 0
//     for k = 1:lc
//       if d(debs(k)-1) > 0
//         if d(fins(k)) < 0
//           imax = [imax round((fins(k)+debs(k))/2)];
//         end
//       else
//         if d(fins(k)) > 0
//           imin = [imin round((fins(k)+debs(k))/2)];
//         end
//       end
//     end
//   end
//
//   if length(imax) > 0
//     indmax = mtlb_sort([indmax imax]);
//   end
//
//   if length(imin) > 0
//     indmin = mtlb_sort([indmin imin]);
//   end
//
// end
// endfunction

//---------------------------------------------------------------------------------------------------

// function ort = emd_io(x,imf)
// // ort = emd_io(x,imf) computes the index of orthogonality
// //
// // inputs : - x    : analyzed signal
// //          - imf  : empirical mode decomposition
//
// n = size(imf,1);
//
// s = 0;
//
// for i = 1:n
//   for j =1:n
//     if i~=j
//       s = s + abs(mtlb_sum(imf(i,:).*conj(imf(j,:)))/mtlb_sum(x.^2));
//     end
//   end
// end
//
// ort = 0.5*s;
// endfunction
//---------------------------------------------------------------------------------------------------

function [x,t,sd,sd2,tol,MODE_COMPLEX,ndirs,display_sifting,sdt,sd2t,r,imf,k,nbit,NbIt,MAXITERATIONS,FIXE,FIXE_H,MAXMODES,INTERP,mask] = init(x,varargin)
    [nargout,nargin]=argn(0);
//x = varargin(1);


// default for stopping
defstop = [0.05,0.5,0.05];

opt_fields = {'t','stop','display','maxiterations','fix','maxmodes','interp','fix_h','mask','ndirs','complex_version'};

defopts.stop = defstop;
defopts.display = 0;
defopts.t = 1:max(size(x));
defopts.maxiterations = 2000;
defopts.fix = 0;
defopts.maxmodes = 0;
defopts.interp = 'spline';
defopts.fix_h = 0;
defopts.mask = 0;
defopts.ndirs = 4;
defopts.complex_version = 2;

opts = defopts;

if size(varargin(1)) > 0
  if typeof(varargin(1)(1))=="st"
    inopts = varargin(1)(1);
    disp(inopts);
  else
  //  error('when using 2 arguments the first one is the analyzed signal X and the second one is a struct object describing the options')
    //inopts = defopts;
    inopts = struct(varargin(1)(1:$));
    disp(inopts);
  end
//elseif nargin > 2
//  try
//    inopts = struct(varargin(2:$));
//  catch
//    error('bad argument syntax')
//  end
  else
  inopts = defopts;
end

if(nargin==1)
  inopts = defopts;
elseif nargin == 0
  error('not enough arguments')
end


fnames = fieldnames(inopts);
for nom = fnames'
  if ~or(convstr(nom) == convstr(opt_fields))
    error(['bad option field name: ',char(nom)])
  end
  if isfield(inopts,nom) // empty values are discarded
     execstr(['opts.'+convstr((nom))+' = inopts.'+(nom)+';'])
  end
end

t = opts.t;
stop = opts.stop;
display_sifting = opts.display;
MAXITERATIONS = opts.maxiterations;
FIXE = opts.fix;
MAXMODES = opts.maxmodes;
INTERP = opts.interp;
FIXE_H = opts.fix_h;
mask = opts.mask;
ndirs = opts.ndirs;
complex_version = opts.complex_version;

if ~(size(x,1)==1 | size(x,2)==1),
  error('X must have only one row or one column')
end

if size(x,1) > 1
  x = x.';
end

if ~(size(t,1)==1 | size(t,2)==1),
  error('option field T must have only one row or one column')
end

if ~isreal(t)
  error('time instants T must be a real vector')
end

if size(t,1) > 1
  t = t';
end

if (length(t)~=length(x))
  error('X and option field T must have the same length')
end

if ~(size(stop,1)==1 | size(stop,2)==1) | length(stop) > 3
  error('option field STOP must have only one row or one column of max three elements')
end

if ~and(abs(x)<%inf)
  error('data elements must be finite')
end

if size(stop,1) > 1
  stop = stop';
end

L = length(stop);
if L < 3
  stop(3)=defstop(3);
end

if L < 2
  stop(2)=defstop(2);
end


if ~type(INTERP)==10 | ~or(convstr(INTERP)==convstr({'linear','cubic','spline'}))
  error('INTERP field must be ''linear'', ''cubic'', ''pchip'' or ''spline''')
end

//special procedure when a masking signal is specified
if or(mask)
  if ~(size(mask,1)==1 | size(mask,2)==1) | length(mask) ~= length(x)
    error('masking signal must have the same dimension as the analyzed signal X')
  end

  if size(mask,1) > 1
    mask = mask.';
  end
  opts.mask = 0;
  imf1 = emd(x+mask,opts);
  imf2 = emd(x-mask,opts);
  if size(imf1,1) ~= size(imf2,1)
    warning('emd: '+['the two sets of IMFs have different sizes: '+string(size(imf1,1))+' and '+string(size(imf2,1))+' IMFs.'])
  end
  S1 = size(imf1,1);
  S2 = size(imf2,1);
  if S1 ~= S2
    if S1 < S2
      tmp = imf1;
      imf1 = imf2;
      imf2 = tmp;
    end
    imf2(max(S1,S2),1) = 0;
  end
  imf = (imf1+imf2)/2;

end


sd = stop(1);
sd2 = stop(2);
tol = stop(3);

lx = length(x);

sdt = sd*ones(1,lx);
sd2t = sd2*ones(1,lx);

if FIXE
  MAXITERATIONS = FIXE;
  if FIXE_H
    error('cannot use both ''FIX'' and ''FIX_H'' modes')
  end
end

if isreal(x)
    MODE_COMPLEX=0
else
  MODE_COMPLEX = complex_version;
end

if MODE_COMPLEX & complex_version ~= 1 & complex_version ~= 2
  error('COMPLEX_VERSION parameter must equal 1 or 2')
end

// number of extrema and zero-crossings in residual
ner = lx;
nzr = lx;

r = x;

if ~or(mask) // if a masking signal is specified "imf" already exists at this stage
  imf = [];
end

k = 1;

// iterations counter for extraction of 1 mode
nbit=0;

// total iterations counter
NbIt=0;
endfunction
