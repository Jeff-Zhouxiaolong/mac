// EMD_SEPARATION.M
//
// P. Flandrin, Mar. 13, 2003 - modified Mar. 2, 2006
//  H. Nahrstaedt - Aug 2010
//
clc;   mode(1);lines(0);
// computes and displays an error measure in the EMD
// separation of two tones
//
// produces Figure 4 in
//
// G. Rilling, P. Flandrin and P. Gon�alv�s
// "On Empirical Mode Decomposition and its algorithms"
// IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
// NSIP-03, Grado (I), June 2003
//halt
mode(-1)
N = 256;// # of data samples
t = 1:N;
tt = N/4:3*N/4;

Nf = 129;// # of tested fequencies
f = linspace(0,.5,Nf);

rapp = 8;// amplitude ratio between modes
a1 = sqrt(rapp);
a2 = 1/a1;

x1 = a1*cos(2*%pi*f'*t);
x2 = a2*cos(2*%pi*f'*t);

se = zeros(Nf,Nf);

for k1 = 5:Nf-1

	for k2 = 2:k1-1;

		y1 = x1(k1,:);
		y2 = x2(k2,:);

		sy1 = sum((y1(tt)).^2);
		sy2 = sum((y2(tt)).^2);
		sy = sum((y1(tt)).^2+(y2(tt)).^2);

		imf = emdc([],y1+y2);
        if size(imf,1) == 1
            imf(2,:) = 0;
        end
		se(k1,k2) = sqrt((sy1*sum((imf(1,tt)-y1(tt)).^2) + sy2*sum((imf(2,tt)-y2(tt)).^2))/(sy1+sy2)/sy);

		//disp([k1 k2 size(imf)]);

	end

end
h=scf(); fig_id=get(gcf(),"figure_id"); clf;
h.color_map = jetcolormap(64);
grayplot(f,f,(se($:-1:1,:)))
//axis('square')
//axis([0 .5 0 .5])
a=gca();a.data_bounds=[0 .5 0 .5];
plot([f(2) f(Nf-1)],[f(Nf-1) f(Nf-1)])
plot([f(Nf-1) f(Nf-1)],[f(2) f(Nf-1)])
plot([f(2) f(Nf-1)],[f(Nf-1) f(2)])
//set(gca,'YTick',[]);set(gca,'XTick',[])
xlabel('f_1 > f_2')
ylabel('f_2')
title('error measure in the EMD separation of two tones');
//colormap(jet);
//colormap(flipud(gray))
colorbar(min(se),max(se))
