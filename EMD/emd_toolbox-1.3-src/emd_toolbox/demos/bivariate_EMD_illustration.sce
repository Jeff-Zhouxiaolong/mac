//bivariate_EMD_illustration.m 
clc;  h=scf(); fig_id=get(gcf(),"figure_id"); clf; mode(1);lines(0);
//illustration of the bivariate EMD extension on a real-world oceanographic signal
//reproduces Fig. 3 in "Bivariate Empirical Mode Decomposition", G. Rilling,
//P. Flandrin, P. Goncalves and J. M. Lilly, IEEE Signal Processing Letters
//
//  H. Nahrstaedt - Aug 2010
//G. Rilling 3/2007 email:  gabriel.rilling@ens-lyon.fr
mode(-1);
demopath = get_absolute_file_path("bivariate_EMD_illustration.sce");
DirectoryStr=demopath+'/data/';
while (isfile([DirectoryStr+'float_position_record.mat'])==0),
 printf('I can''t find %s\n', [DirectoryStr 'float_position_record.mat']);
 DirectoryStr=input('name of the directory where float_position_record.mat is : ','s');
end;
loadmatfile(DirectoryStr+'float_position_record.mat');

//loadmatfile('data/float_position_record.mat');

[imf,nb] = cemdc2_fix([],x,10,[],32);
n = size(imf,1);

figtitle1 = 'Float position record';
plot(real(x),imag(x));
xlabel('Displacement East (km) --- Real part')
ylabel('Displacement North (km) --- Imaginary part')
title(figtitle1)
//axis equal;
//set(gca,'Ylim',[-250,300])

figtitle2 = 'Bivariate Empirical Mode Decomposition of Float signal';
scf()
subplot(n+1,1,1)
plot(real(x))
//hold on
plot(imag(x),'k--')
//axis tight
ylabel('signal')
title(figtitle2)
//set(gca,'XTickLabel',{})

m = min([real(imf(1:$-1,:));imag(imf(1:$-1,:))]);
M = max([real(imf(1:$-1,:));imag(imf(1:$-1,:))]);
for k = 1:n
  subplot(n+1,1,k+1)
  plot(real(imf(k,:)))
  //hold on
  plot(imag(imf(k,:)),'k--')
  if k<n
    a=gca();a.data_bounds=([1,length(x),m,M]);
  end
  ylabel(['d_'+string(k)])
  //if k<n
  //  set(gca,'XTickLabel',{})
  //end
end
ylabel('res.')
xlabel('Time (days)')
//axis tight
