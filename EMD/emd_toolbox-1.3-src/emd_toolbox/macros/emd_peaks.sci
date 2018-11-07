function [pks_max,locs_max] =emd_peaks(X)
//  finds peaks (local maxima)
// Calling Sequence
// [pks_max,locs_max] =emd_peaks(X)
// Parameters
// inputs :
//          - X : analyzed signal
// outputs :
//           - pks_max  = peaks
//           - locs_max = indices of peaks
//
// See also
//  boundary_conditions_emd
// Authors
// H. Nahrstaedt - Aug 2010

nmax = 0;                  // counter for max peaks
L = length(X);
j = 0;
pks_max  = zeros(1,L);
locs_max = zeros(1,L);
for j=2:L-1
    if(and((X(j) > [X(j-1) X(j+1)])))
        nmax = nmax+1;
        pks_max(nmax)  = X(j);
        locs_max(nmax) = j;
    end
end
if nmax~=0
    pks_max  = pks_max(1:nmax);
    locs_max = locs_max(1:nmax);
else
    pks_max  = [];
    locs_max = [];
end
endfunction
