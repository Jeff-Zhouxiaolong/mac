function [indmin, indmax] = emd_local_peaks(x)
//  finds local extrema
// Calling Sequence
// [indmin, indmax] = emd_local_peaks(x)
// Parameters
// inputs :
//          - x : analyzed signal
// outputs :
//           - indmin = indices of minima
//           - indmax = indices of maxima
// Examples
// x=sin(2*%pi*10*linspace(0,1,100));
// [indmin, indmax] = emd_local_peaks(x);
// scf(); plot(x);
// plot(indmin,x(indmin),'or');
// plot(indmax,x(indmax),'oc');
// See also
//  peaks
//  emd_zero_crossings
//  boundary_conditions_emd
// Authors
// H. Nahrstaedt - Aug 2010

// if(and(x < 1e-5))
//     x=zeros(1,length(x));
// end
if ~isempty(x)
m = length(x);
// Calculates the extrema of the projected signal
// Difference between subsequent elements:
dy = diff(x); a = find(dy~=0);
lm = find(diff(a)~=1) + 1;
if lm>1 & length(a)>=lm
  d = a(lm) - a(lm-1);
  a(lm) = a(lm) - floor(d/2);
  a($+1) = m;
else
  a($+1)=m;
end
ya  = x(a);
else
ya=[];
end;
if(length(ya) > 1)
    // Maxima
    [pks_max,loc_max]=emd_peaks(ya);
    // Minima
    [pks_min,loc_min]=emd_peaks(-1*ya);

    if(~isempty(pks_min))
        indmin = a(loc_min);
    else
        indmin = %nan;
    end
    if(~isempty(pks_max))
        indmax = a(loc_max);
    else
        indmax = %nan;
    end
else
    indmin=%nan;
    indmax=%nan;
end
endfunction
