function [f2c, c2f] = emd_reconstruct(imf);
 //partial reconstructions (fine to coarse and coarse to fine) for real and complex data.
// Calling Sequence
//   [f2c, c2f] = emd_reconstruct(imf)
// Parameters
// inputs :   
//            - imf : intrinsic mode functions (last line = residual)
// outputs : 
//            - f2c : fine to coarse reconstruction
//            - c2f : coarse to fine reconstruction
// Examples 
//     s = rand(1,512,'normal');
//     imf = emd(s);
//     [f2c, c2f] = emd_reconstruct(imf);
// See also
//  cemd_visu
//  emd_visu
// Authors
// H. Nahrstaedt - Aug  2010
// P. Flandrin, Mar. 13, 2003
// G. Rilling, last modification 3.2006 gabriel.rilling@ens-lyon.fr

s = size(imf);
k = s(1);

f2c = [];
f2c(1,:) = imf(1,:);

c2f = [];
c2f(1,:) = imf(k,:);

for j = 2:k
  f2c(j,:) = f2c(j-1,:) + imf(j,:);
  c2f(j,:) = c2f(j-1,:) + imf(k+1-j,:);
end

endfunction