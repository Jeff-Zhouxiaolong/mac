// EMD_SAMPLING.M
//
// P. Flandrin, Mar. 13, 2003 - modified Mar. 2, 2006
//  H. Nahrstaedt - Aug 2010
//
clc;  scf(); fig_id=get(gcf(),"figure_id"); clf; mode(1);lines(0);
// computes and plots an error measure in the EMD
// estimation of a single tone
//
// produces Figure 3 in
//
// G. Rilling, P. Flandrin and P. Gon�alves
// "On Empirical Mode Decomposition and its algorithms"
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003

mode(-1)
N = 256;// # of data samples
t = 1:N;
tt = fix(N/4):fix(3*N/4);

Nf = 257;// # of tested fequencies
f = logspace(-log10(2*Nf),-log10(2),Nf);

x = cos(2*%pi*f'*t);

se = zeros(1,Nf);

kmin = 65;

for k = kmin:Nf-1

	y = x(k,:);

	sy = sum((y(tt)).^2);

	imf = emdc([],y);
	se(k) = sqrt(sum((imf(1,tt)-y(tt)).^2)/sy);

	//disp([k f(k) size(imf,1)]);

end
tmp = se(kmin:Nf-1);
for k=1:length(tmp)
   if tmp(k)==0
     tmp(k)=-%inf;
   else
     tmp(k)=log2(tmp(k));
   end
end

plot(log2(f(kmin:Nf-1)),max(tmp,-60),'o-')
//axis([-8 -1 -16 0])
a=gca();a.data_bounds=[-8 -1 -16 0];
xlabel('log_{2}(frequency)')
ylabel('log_{2}(error)')

plot([-8 -1],[-14 0],'r')
xgrid
