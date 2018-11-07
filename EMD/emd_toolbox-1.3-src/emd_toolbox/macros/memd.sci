function q = memd(x, varargin)
//Multivariate Empirical Mode Decomposition algorithm for signals cotaining 3-16 channels.
// Calling Sequence
//imf = memd(X)
//imf = memd(X,num_directions)
//imf = memd(X,num_directions,'stopping criteria'
// imf = memd(X, num_directions, 'stop', stop_vec)
//imf = memd(X, num_directions, 'fix_h', n_iter)
// Description
// function NEMD applies the "Multivariate Empirical Mode Decomposition" algorithm (Rehman and Mandic, Proc. Roy. Soc A, 2010)
// to multivariate inputs.
//   We have verified this code by simulations for signals cotaining 3-16 channels.
//
// imf = NEMD(X)
//   returns a 3D matrix 'imf(N,M,L)' containing M multivariate IMFs, one IMF per column, computed by applying
//   the multivariate EMD algorithm on the N-variate signal (time-series) X of length L.
//    - For instance, imf_k = IMF(k,:,:) returns the k-th component (1 <= k <= N) for all of the N-variate IMFs.
//
//   For example,  for hexavariate inputs (N=6), we obtain a 3D matrix IMF(6, M, L)
//   where M is the number of IMFs extracted, and L is the data length.
//
// imf = NEMD(X,num_directions)
//   where integer variable num_directions (>= 1) specifies the total number of projections of the signal
//     - As a rule of thumb, the minimum value of num_directions should be twice the number of data channels,
//     - for instance, num_directions = 6  for a 3-variate signal and num_directions= 16 for an 8-variate signal
//   The default number of directions is chosen to be 128 - to extract meaningful IMFs, the number of directions
//   should be considerably greater than the dimensionality of the signals
//
// imf = NEMD(X,num_directions,'stopping criteria')
//   uses the optional parameter 'stopping criteria' to control the sifting process.
//    The available options are
//      -  'stop' which uses the standard stopping criterion specified in [2]
//      -  'fix_h' which uses the modified version of the stopping criteria specified in [3]
//    The default value for the 'stopping criteria' is 'stop'.
//
//  The settings  num_directions=128 and 'stopping criteria' = 'stop' are defaults.
//     Thus imf = NEMD(X) = NEMD(X,128) = NEMD(X,128,'stop') = NEMD(X,[],'stop'),
//
// imf = NEMD(X, num_directions, 'stop', stop_vec)
//   computes the IMFs based on the standard stopping criterion whose parameters are given in the 'stop_vec'
//     - stop_vec has three elements specifying the threshold and tolerance values used, see [2].
//     - the default value for the stopping vector is   step_vec = [0.075 0.75 0.075].
//     - the option 'stop_vec' is only valid if the parameter 'stopping criteria' is set to 'stop'.
//
// imf = NEMD(X, num_directions, 'fix_h', n_iter)
//   computes the IMFs with n_iter (integer variable) specifying the number of consecutive iterations when
//   the number of extrema and the number of zero crossings differ at most by one [3].
//     - the default value for the parameter n_iter is set to  n_iter = 5.
//     - the option n_iter is only valid if the parameter  'stopping criteria' = 'fix_h'
//
//
// This code allows to process multivaraite signals having 3-16 channels, using the multivariate EMD algorithm [1].
//   - to perform EMD on more than 16 channels, modify the variable 'Max_channels' on line 536 in the code accordingly.
//   - to process 1- and 2-dimensional (univariate and bivariate) data using EMD, we recommend the toolbox from
//                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html
//
// Acknowledgment: Part of this code is based on the bivariate EMD code, publicly available from
//                 http://perso.ens-lyon.fr/patrick.flandrin/emd.html
//
//
// Bibliography
//
// [1]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010
//
// [2]  G. Rilling, P. Flandrin and P. Gonçalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
//      Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003
//
// [3]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
//      Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003
//
// Examples
// Case 1:
//
// inp = rand(1000,3,'normal');
// imf = memd(inp);
// imf_x = matrix(imf(1,:,:).entries,size(imf,2),size(imf,3));  // imfs corresponding to 1st component
// imf_y = matrix(imf(2,:,:).entries,size(imf,2),size(imf,3));  // imfs corresponding to 2nd component
// imf_z = matrix(imf(3,:,:).entries,size(imf,2),size(imf,3));  // imfs corresponding to 3rd component
//
// Case 2:
//
// loadmatfile(emd_getpath()+"demos/data/syn_hex_inp.mat");
// imf = memd(s6,256,'stop',[0.05 0.5 0.05])
// Authors
// H. Nahrstaedt - Aug 2010
// Copyright: Naveed ur Rehman and Danilo P. Mandic, Oct-2009

[nargout,nargin]=argn(0);
global N N_dim;
[x, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stop_crit, stp_cnt] = set_value(x, nargin, varargin);

r=x; n_imf=1;
while ~stop_emd(r, seq, ndir)
    // current mode
    m = r;
    // mode at previous iteration
    mp = m;

    // computation of mean and stopping criterion
    if((stop_crit=='stop'))
        [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
    else
        counter=0;
        [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
    end

    // In case the current mode is so small that machine precision can cause
    // spurious extrema to appear
    if (max(abs(m))) < (1e-10)*(max(abs(x)))
        if ~stop_sift
            warning('emd:warning','forced stop of EMD : too small amplitude')
        else
            disp('forced stop of EMD : too small amplitude')
        end
        break
    end

    // sifting loop
    while ~stop_sift & nbit<MAXITERATIONS
        //sifting
        m = m - env_mean;
        // computation of mean and stopping criterion
        if((stop_crit=='stop'))
            [stop_sift,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir);
        else
            [stop_sift,env_mean,counter] = stop_sifting_fix(m,t,seq,ndir,stp_cnt,counter);
        end
        mp = m;
        nbit=nbit+1;

        if(nbit==(MAXITERATIONS-1) &  nbit > 100)
                warning('emd:warning',['forced stop of sifting : too many iterations']);
        end
    end
    for n_dim=1:N_dim
        q(n_dim,n_imf,:)=m(:,n_dim);
    end
    n_imf = n_imf+1;
    r = r - m;
    nbit = 0;
end
// Stores the residue
for n_dim=1:N_dim
    q(n_dim,n_imf,:)=r(:,n_dim);
end

//sprintf('Elapsed time: //f\n',toc);
endfunction

//---------------------------------------------------------------------------------------------------
function stp = stop_emd(r, seq, ndir)
global N N_dim;
for it=1:ndir
    if (N_dim~=3) // Multivariate signal (for N_dim ~=3) with hammersley sequence
        // Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b=2*seq(1:$,it)-1;

        // Find angles corresponding to the normalised sequence
        for n=1:N_dim-1
            tht(n)=atan(norm(b(n+1:N_dim),2),b(n));
        end
        // Find coordinates of unit direction vectors on n-sphere
        for n=1:N_dim
            if n == N_dim
                dir_vec(n)=1;
            else
                dir_vec(n)=cos(tht(n));
            end
            for i=1:n-1
                dir_vec(n)=dir_vec(n)*sin(tht(i));
            end
        end
    else // Trivariate signal with hammersley sequence
        // Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt(find(tt>1))=1;
        tt(find(tt<-1))=-1;

        // Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*%pi;
        st = sqrt(1.0-tt*tt);

        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end
    // Projection of input signal on nth (out of total ndir) direction vectors
    for iter=1:N
        y(iter)=dot(dir_vec, r(iter,:));
    end
    // Calculates the extrema of the projected signal
    [indmin, indmax] = emd_local_peaks(y);


    ner(it) = length(indmin) + length(indmax);
end

// Stops if the all projected signals have less than 3 extrema
stp = and(ner < 3);
endfunction

//---------------------------------------------------------------------------------------------------
// computes the mean of the envelopes and the mode amplitude estimate
function [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir) //new
global N N_dim;
NBSYM = 2;
count=0;

env_mean=zeros(length(t),N_dim);
env_min=zeros(length(t),N_dim);
env_max=zeros(length(t),N_dim);
amp = zeros(length(t),1);

for it=1:ndir
    if (N_dim ~=3) // Multivariate signal (for N_dim ~=3) with hammersley sequence
        // Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b=2*seq(1:$,it)-1;
        // Find angles corresponding to the normalised sequence
        for n=1:N_dim-1
            tht(n)=atan(norm(b(n+1:N_dim),2),b(n));
        end
        // Find coordinates of unit direction vectors on n-sphere
        for n=1:N_dim
            if n == N_dim
                dir_vec(n)=1;
            else
                dir_vec(n)=cos(tht(n));
            end
            for i=1:n-1
                dir_vec(n)=dir_vec(n)*sin(tht(i));
            end
        end
    else // Trivariate signal with hammersley sequence
        // Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt(find(tt>1))=1;
        tt(find(tt<-1))=-1;

        // Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*%pi;
        st = sqrt(1.0-tt*tt);

        dir_vec(1)=st * cos(phirad);
        dir_vec(2)=st * sin(phirad);
        dir_vec(3)=tt;
    end

    // Projection of input signal on nth (out of total ndir) direction vectors
    for iter=1:N
        y(iter)=dot(dir_vec, m(iter,:));
    end

    // Calculates the extrema of the projected signal
    [indmin, indmax] = emd_local_peaks(y);


    nem(it) = length(indmin) + length(indmax);

    indzer = emd_zero_crossings(y);
    nzm(it) = length(indzer);

    [tmin,tmax,zmin,zmax,emode] = boundary_conditions_emd(indmin,indmax,t,y,m,NBSYM);

    // Calculate multidimensional envelopes using spline interpolation
    // Only done if number of extrema of the projected signal exceed 3
    if(emode)
        for n=1:N_dim
            env_min(:,n)=interp1(tmin,zmin(:,n),t,'spline')';
            env_max(:,n)=interp1(tmax,zmax(:,n),t,'spline')';
        end
        for iter=1:length(amp)
            amp(iter) = amp(iter) + norm(env_max(iter,:)-env_min(iter,:))/2;
        end
        env_mean = env_mean + (env_max+env_min)/2;
    else // if the projected signal has inadequate extrema
        count=count+1;
    end
end
if(ndir>count)
    env_mean = env_mean/(ndir-count);
    amp = amp/(ndir-count);
else
    env_mean = zeros(N,N_dim);
    amp = zeros(N,1);
    nem = zeros(1,ndir);
end
endfunction

//-------------------------------------------------------------------------------
// Stopping criterion
function [stp,env_mean] = stop_sifting(m,t,sd,sd2,tol,seq,ndir)
global N N_dim;
try
    [env_mean,nem,nzm,amp] = envelope_mean(m,t,seq,ndir);
    if(amp)
        for it=1:length(amp)
            sx(it) = norm(env_mean(it,:))./amp(it);
        end
    else
        for it=1:length(amp)
            sx(it) = norm(env_mean(it,:));
        end
    end
    stp = ~((mean(bool2s(sx > sd)) > tol | or(sx > sd2)) & or(nem > 2));
catch
    env_mean = zeros(N,N_dim);
    stp = 1;
end
endfunction

function [stp,env_mean,counter]= stop_sifting_fix(m,t,seq,ndir,stp_count,counter)
global N N_dim;
try
[env_mean,nem,nzm] = envelope_mean(m,t,seq,ndir);
if (and(abs(nzm-nem)>1))
    stp = 0;
    counter = 0;
else
    counter = counter+1;
    stp = (counter >= stp_count);
end
catch
  env_mean = zeros(N,N_dim);
  stp = 1;
end
endfunction

//---------------------------------------------------------------------------------------
// defines new extrema points to extend the interpolations at the edges of the
// signal (mainly mirror symmetry)
// function [tmin,tmax,zmin,zmax,emode] = boundary_conditions(indmin,indmax,t,x,z,nbsym)
//
// lx = length(x);
// if (length(indmin) + length(indmax) < 3)
//     emode = 0;
//     tmin=%nan;tmax=%nan;zmin=%nan;zmax=%nan;
//     return
// else
//     emode=1; //the projected signal has inadequate extrema
// end
//     // boundary conditions for interpolations :
//     if indmax(1) < indmin(1)
//     	if x(1) > x(indmin(1))
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
//     end
//     tlmin = 2*t(lsym)-t(lmin);
//     tlmax = 2*t(lsym)-t(lmax);
//     trmin = 2*t(rsym)-t(rmin);
//     trmax = 2*t(rsym)-t(rmax);
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
//     if trmin($) < t(lx) | trmax($) < t(lx)
//         if rsym == indmax($)
//             rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
//         else
//             rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
//         end
//         if rsym == lx
//             error('bug')
//         end
//         rsym = lx;
//         trmin = 2*t(rsym)-t(rmin);
//         trmax = 2*t(rsym)-t(rmax);
//     end
// 	zlmax =z(lmax,:);
// 	zlmin =z(lmin,:);
// 	zrmax =z(rmax,:);
// 	zrmin =z(rmin,:);
//
// 	tmin = [tlmin t(indmin) trmin];
// 	tmax = [tlmax t(indmax) trmax];
// 	zmin = [zlmin; z(indmin,:); zrmin];
// 	zmax = [zlmax; z(indmax,:); zrmax];
// endfunction
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

//
// function [indmin, indmax] = emd_local_peaks(x)
// if(and(x < 1e-5))
//     x=zeros(1,length(x));
// end
// m = length(x);
// // Calculates the extrema of the projected signal
// // Difference between subsequent elements:
// dy = diff(x); a = find(dy~=0);
// lm = find(diff(a)~=1) + 1;
// if lm>1
//   d = a(lm) - a(lm-1);
//   a(lm) = a(lm) - floor(d/2);
//   a($+1) = m;
// else
//   a($+1)=m;
// end
// ya  = x(a);
//
// if(length(ya) > 1)
//     // Maxima
//     [pks_max,loc_max]=emd_peaks(ya);
//     // Minima
//     [pks_min,loc_min]=emd_peaks(-1*ya);
//
//     if(~isempty(pks_min))
//         indmin = a(loc_min);
//     else
//         indmin = %nan;
//     end
//     if(~isempty(pks_max))
//         indmax = a(loc_max);
//     else
//         indmax = %nan;
//     end
// else
//     indmin=%nan;
//     indmax=%nan;
// end
// endfunction
//
// function [pks_max,locs_max] =emd_peaks(X)
// nmax = 0;                  // counter for max peaks
// L = length(X);
// j = 0;
// pks_max  = zeros(1,L);
// locs_max = zeros(1,L);
// for j=2:L-1
//     if(and((X(j) > [X(j-1) X(j+1)])))
//         nmax = nmax+1;
//         pks_max(nmax)  = X(j);
//         locs_max(nmax) = j;
//     end
// end
// if nmax~=0
//     pks_max  = pks_max(1:nmax);
//     locs_max = locs_max(1:nmax);
// else
//     pks_max  = [];
//     locs_max = [];
// end
// endfunction
//
// function indzer = emd_zero_crossings(x)
//   indzer = find(x(1:$-1).*x(2:$)<0);
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
// endfunction

function seq = hamm(n,base)
  seq = zeros(1,n);
    if ( 1 < base )
      seed = 1:1:n;
      base_inv = inv(base);
      while ( or ( seed ~= 0 ) )
        digit = pmodulo (seed(1:n), base);
        seq = seq + digit * base_inv;
        base_inv = base_inv / base;
        seed = floor (seed / base );
      end
    else
      temp = 1:1:n;
      seq = (pmodulo(temp,(-base + 1 ))+0.5)/(-base);
    end
endfunction

function [q, seq, t, ndir, N_dim, N, sd, sd2, tol, nbit, MAXITERATIONS, stp_crit, stp_cnt] = set_value(q, narg, varargin)

//error(nargchk(1,4,narg));
if (narg==0)
  error("At least one parameter needed!");
end



if  size(varargin(1)) > 0 & type(varargin(1))==15
varargin = varargin(1);
end;
//disp(varargin)

ndir = [];
stp_crit = [];
stp_vec = [];
stp_cnt  = [];
MAXITERATIONS  = [];
sd=[];
sd2=[];
tol=[];
prm= [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149];

// Changes the input vector to double vector
q = double(q);

// Specifies maximum number of channels that can be processed by the code
// Its maximum possible value is 32.
Max_channels = 16;

if(narg==2)
    ndir=varargin(1);
end

if(narg==3)
    if(~isempty(varargin(1)))
        ndir=varargin(1);
    else
        ndir=128;
    end
    stp_crit=varargin(2);
end

if(narg==4 & (varargin(2)=='fix_h'))
    if(isempty(varargin(1)))
        ndir=128;
        stp_crit=varargin(2);
        stp_cnt  = varargin(3);
    else
        ndir=varargin(1);
        stp_crit=varargin(2);
        stp_cnt  = varargin(3);
    end
elseif (narg==4 & (varargin(2)=='stop'))
    if(isempty(varargin(1)))
        ndir=128;
        stp_crit=varargin(2);
        stp_vec=varargin(3);
    else
        ndir=varargin(1);
        stp_crit=varargin(2);
        stp_vec=varargin(3);
    end
elseif (narg==4 & ((varargin(2)~='fix_h') & (varargin(2)~='stop')))
    //Nmsgid = generatemsgid('invalid stop_criteria');
    error('stop_criteria should be either fix_h or stop');
end

//////////////////////////// Rescale input signal if required
if (or(size(q)) == 0)
    //datamsgid = generatemsgid('emptyDataSet');
    error('Data set cannot be empty.');
end
if size(q,1) < size(q,2)
   q=q';
end

//////////////////////// Dimension of input signal
N_dim = size(q,2);
if(N_dim < 3 | N_dim > Max_channels)
    error('Function only processes the signal having 3 and 16 channels.');
end

//////////////////////// Length of input signal
N = size(q,1);

////////////////////////// Check validity of Input parameters
if ~isempty(ndir) & (~or(type(ndir)==[1 5 8]) | ~(sum(length(ndir))==1) | or(ndir-fix(ndir./1).*1) | (ndir < 6))
    //Nmsgid = generatemsgid('invalid num_dir');
    error('num_dir should be an integer greater than or equal to 6.');
end

if ~isempty(stp_crit) & (~ischar(stp_crit) | ((stp_crit~='fix_h')&(stp_crit~='stop')))
   // Nmsgid = generatemsgid('invalid stop_criteria');
    error('stop_criteria should be either fix_h or stop');
end

if ~isempty(stp_vec) & (~or(type(stp_vec)==[1 5 8]) | length(stp_vec)~=3 | ~(stp_crit=='stop'))
   // Nmsgid = generatemsgid('invalid stop_vector');
    error('stop_vector should be an array with three elements e.g. default is [0.075 0.75 0.075] ');
end

if ~isempty(stp_cnt) & (~or(type(stp_cnt)==[1 5 8]) | ~(sum(length(stp_cnt))==1) | or(stp_cnt-fix(stp_cnt./1).*1) | (stp_cnt < 0) | ~(stp_crit=='fix_h'))
   // Nmsgid = generatemsgid('invalid stop_count');
    error('stop_count should be a nonnegative integer');
end

if (isempty(ndir))
    ndir=128; // default
end

if (isempty(stp_crit))
    stp_crit='stop'; // default
end

if (isempty(stp_vec))
    stp_vec=[0.075,0.75,0.075]; // default
end

if (isempty(stp_cnt))
    stp_cnt=5; // default
end

if((stp_crit=='stop'))
    sd = stp_vec(1);
    sd2 = stp_vec(2);
    tol = stp_vec(3);
end

////////////////////////// Initializations for Hammersley function
base(1) = -ndir;

//////////////////////////// Find the pointset for the given input signal
    if(N_dim==3)
        base(2) = 2;
        for it=1:N_dim-1
            seq(it,:) = hamm(ndir,base(it));
        end
        //seq = i4_to_hammersley_sequence (2, ndir, base );
    else
        for iter = 2 : N_dim
            base(iter) = prm(iter-1);
        end

        for it=1:N_dim
            seq(it,:) = hamm(ndir,base(it));
        end
    end

//////////////////////// Define t
t=1:N;

// Counter
nbit=0;
MAXITERATIONS=1000; // default

// disp([ndir,stp_crit,stp_vec,stp_cnt,MAXITERATIONS,sd,sd2,tol]);
//disp(sprintf('ndir: //d\nstop_criteria: //s\nstop_vec(in case of "stop" criteria): [//d //d //d]\nstop_count(in case of "fix" criteria)://d', ndir, stp_crit, sd, sd2, tol,stp_cnt));
// tic
endfunction

function z = dot(x, y, dim)
 [nargout,nargin]=argn(0);
  if (nargin ~= 2 & nargin ~= 3)
  error("wrong number of parameters!");
  end

  if (nargin < 3)
    if  size(x,1)==1 | size(x,2)==1
      x = x(:);
    end
    if size(y,1)==1 | size(y,2)==1
      y = y(:);
    end
    if (abs(max(size(x)-size(y)))>0 )
      error ("dot: sizes of arguments must match");
    end
    z = sum(conj (x) .* y);
  else
    if (abs(max(size(x)-size(y)))>0 )
      error ("dot: sizes of arguments must match");
    end
    z = sum(conj (x) .* y, dim);
  end

endfunction
