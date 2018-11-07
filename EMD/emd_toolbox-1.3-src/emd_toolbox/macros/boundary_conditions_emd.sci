function [tmin,tmax,zmin,zmax,emode] = boundary_conditions_emd(indmin,indmax,t,x,z,nbsym)
//  extends an extrema set to limit edge effects on the interpolations
// Calling Sequence
// [TMIN,TMAX,ZMIN,ZMAX,EMODE] = boundary_conditions_emd(INDMIN,INDMAX,T,X,Z,NBSYM)
// [TMIN,TMAX,ZMIN,ZMAX,EMODE] = boundary_conditions_emd(INDMIN,INDMAX,[],X,Z,NBSYM)
// Parameters
// inputs:
//   - INDMIN, INDMAX: indices of minima and maxima in the real signal X
//   - T: sampling times. If emtpy, then t is set to t = 1:length(x);
//   - X: real signal in which INDMIN and INDMAX are the indices of extrema
//   - Z: signal which values are interpolated in the final envelope
//   - NBSYM: number of points added to each end
// outputs:
//   - TMIN, TMAX: extended sampling times
//   - ZMIN, ZMAX: extended "extrema" set
//  -  EMODE: 0 means the signal has not enough extrema. 1 means x has enough extrema.
// Description
// defines new extrema points to extend the interpolations at the edges of the
// signal (mainly mirror symmetry)
//
//   - for a real signal X:
//
//     [TMIN,TMAX,ZMIN,ZMAX] = boundary_conditions_emd(INDMIN,INDMAX,T,X,X,NBSYM)
//
//   - for a complex signal Z and a direction PHI:
//
//     X = exp(-i*PHI)*Z;
//
//     [TMIN,TMAX,ZMIN,ZMAX] = boundary_conditions_emd(INDMIN,INDMAX,T,X,Z,NBSYM)
//
// TODO: it has to be noted that this function was originally written for the
// classical EMD and adapted to the bivariate case without a proper study of its
// effects. The edge effects problem for the bivariate EMD has not been studied yet.
// See also
//  emd_local_peaks
// Authors
// H. Nahrstaedt
// G. Rilling, last modification 3.2007 gabriel.rilling@ens-lyon.fr


 [nargout,nargin]=argn(0);
  if (nargin<6)
     error("At least 6 parameter needed!");
  end
	lx = length(x);
	if (length(indmin) + length(indmax) < 3)
		//error('not enough extrema')
	      emode = 0;
	      tmin=%nan;tmax=%nan;zmin=%nan;zmax=%nan;
	      return
	else
	    emode=1; //the projected signal has inadequate extrema
	end
      if isempty(t)
	t = 1:length(x);
      end
    // boundary conditions for interpolations :

	if indmax(1) < indmin(1)
    	if x(1) > x(indmin(1))
			lmax = moc_fliplr(indmax(2:min(length(indmax),nbsym+1)));
			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
			lsym = indmax(1);
		else
			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
			lmin = [moc_fliplr(indmin(1:min(length(indmin),nbsym-1))),1];
			lsym = 1;
		end
	else

		if x(1) < x(indmax(1))
			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
			lmin = moc_fliplr(indmin(2:min(length(indmin),nbsym+1)));
			lsym = indmin(1);
		else
			lmax = [moc_fliplr(indmax(1:min(length(indmax),nbsym-1))),1];
			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
			lsym = 1;
		end
	end

	if indmax($) < indmin($)
		if x($) < x(indmax($))
			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym,1):$-1));
			rsym = indmin($);
		else
			rmax = [lx,moc_fliplr(indmax(max(length(indmax)-nbsym+2,1):$))];
			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
			rsym = lx;
		end
	else
		if x($) > x(indmin($))
			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym,1):$-1));
			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
			rsym = indmax($);
		else
			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
			rmin = [lx,moc_fliplr(indmin(max(length(indmin)-nbsym+2,1):$))];
			rsym = lx;
		end
	end

	tlmin = 2*t(lsym)-t(lmin);
	tlmax = 2*t(lsym)-t(lmax);
	trmin = 2*t(rsym)-t(rmin);
	trmax = 2*t(rsym)-t(rmax);

	// in case symmetrized parts do not extend enough
	if tlmin(1) > t(1) | tlmax(1) > t(1)
		if lsym == indmax(1)
			lmax = moc_fliplr(indmax(1:min(length(indmax),nbsym)));
		else
			lmin = moc_fliplr(indmin(1:min(length(indmin),nbsym)));
		end
		if lsym == 1
			error('bug')
		end
		lsym = 1;
		tlmin = 2*t(lsym)-t(lmin);
		tlmax = 2*t(lsym)-t(lmax);
	end

	if trmin($) < t(lx) | trmax($) < t(lx)
		if rsym == indmax($)
			rmax = moc_fliplr(indmax(max(length(indmax)-nbsym+1,1):$));
		else
			rmin = moc_fliplr(indmin(max(length(indmin)-nbsym+1,1):$));
		end
		if rsym == lx
			error('bug')
		end
		rsym = lx;
		trmin = 2*t(rsym)-t(rmin);
		trmax = 2*t(rsym)-t(rmax);
	end
      if (size(z,1)>1)
	zlmax =z(lmax,:);
	zlmin =z(lmin,:);
	zrmax =z(rmax,:);
	zrmin =z(rmin,:);

	tmin = [tlmin t(indmin) trmin];
	tmax = [tlmax t(indmax) trmax];
	zmin = [zlmin z(indmin,:) zrmin];
	zmax = [zlmax z(indmax,:) zrmax];
     else
	zlmax =z(lmax);
	zlmin =z(lmin);
	zrmax =z(rmax);
	zrmin =z(rmin);

	tmin = [tlmin t(indmin) trmin];
	tmax = [tlmax t(indmax) trmax];
	zmin = [zlmin z(indmin) zrmin];
	zmax = [zlmax z(indmax) zrmax];
     end

endfunction
