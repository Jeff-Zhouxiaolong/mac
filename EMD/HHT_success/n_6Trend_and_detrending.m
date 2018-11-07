load gsta.dat;

year=gsta(:,1);

inData=gsta(:,2);

rslt=eemd(inData,0,1);

plot(year, rslt(:,1));

hold on;

plot(year, sum(rslt(:,7:8),2),'r-');

plot(year, sum(rslt(:,6:8),2),'g-');

plot(year, sum(rslt(:,5:8),2),'m-');

title('Trends of different timescales');

ylabel('Kelvin');

xlabel('year');

grid;

