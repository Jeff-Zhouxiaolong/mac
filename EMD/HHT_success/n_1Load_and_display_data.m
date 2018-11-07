load gsta.dat;

plot(gsta(:,1),gsta(:,2));

axis([1850 2010 -0.6 0.6]);

title('the annual mean global surface temperature anomaly')

xlabel('year')

ylabel('Kelvin')

