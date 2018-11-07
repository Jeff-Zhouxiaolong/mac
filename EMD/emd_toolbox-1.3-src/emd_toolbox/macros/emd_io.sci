function ort = emd_io(x,imf)
//  computes the index of orthogonality
// Calling Sequence
// ORT = emd_io(X,IMF)
// Parameters
// inputs :
//           - X    : analyzed signal
//          - IMF  : empirical mode decomposition
// Authors
// H. Nahrstaedt - Aug 2010
// G. Rilling, last modification: 3.2007

n = size(imf,1);

s = 0;

for i = 1:n
  for j =1:n
    if i~=j
      s = s + abs(mtlb_sum(imf(i,:).*conj(imf(j,:)))/mtlb_sum(x.^2));
    end
  end
end

ort = 0.5*s;
endfunction
