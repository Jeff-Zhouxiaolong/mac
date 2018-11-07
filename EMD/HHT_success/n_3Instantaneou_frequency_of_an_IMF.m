omega_m3=ifndq(rslt(:,4),1);

subplot(2,1,1);

plot(year,rslt(:,4));

axis([1850 2010 -0.12 0.12]); 

title('IMF C3');

ylabel('Kelvin');

grid;

subplot(2,1,2);

plot(year, omega_m3/2/pi,'r-');

grid;

xlabel('year');

ylabel('cycle/year');

title('instantaneous frequency');

axis([1850 2010 0 0.12]);

