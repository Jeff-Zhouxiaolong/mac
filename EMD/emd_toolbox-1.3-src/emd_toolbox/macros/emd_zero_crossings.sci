function indzer = emd_zero_crossings(x)
//  finds zero-crossings
// Calling Sequence
// indzer = emd_zero_crossings(x)
// Parameters
// inputs :
//          - x : analyzed signal
// outputs :
//           - indzer: = indices of zero-crossings
// Examples
// x=sin(2*%pi*10*linspace(0,1,100));
// [indzer] = emd_zero_crossings(x);
// scf(); plot2d2(x);
// plot(indzer,x(indzer),'or');
// See also
//  boundary_conditions_emd
// Authors
// H. Nahrstaedt - Aug 2010

  indzer = find(x(1:$-1).*x(2:$)<0);
  if or(x == 0)
    iz = find( x==0 );
    indz = [];
    if or(diff(iz)==1)
      zer = x == 0;
      dz = diff([0 zer 0]);
      debz = find(dz == 1);
      finz = find(dz == -1)-1;
      indz = round((debz+finz)/2);
    else
      indz = iz;
    end
    indzer = mtlb_sort([indzer indz]);
  end
endfunction
