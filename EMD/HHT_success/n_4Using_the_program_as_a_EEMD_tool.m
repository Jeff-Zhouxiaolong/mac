load gsta.dat;

year=gsta(:,1);

inData=gsta(:,2);

rslt=eemd(inData,0,1);

rslt=eemd(inData,0.2,100);

t(1)=1850;

t(2)=2010;

y1(1)=0;

y1(2)=0;

y2(1)=-0.3;

y2(2)=-0.3;

y3(1)=-0.6;

y3(2)=-0.6;

y4(1)=-0.9;

y4(2)=-0.9;

y5(1)=-1.2;

y5(2)=-1.2;

y6(1)=-1.6;

y6(2)=-1.6;

plot(t,y1,'k-');

hold on;

plot(t,y2,'k-');

plot(t,y3,'k-');

plot(t,y4,'k-');

plot(t,y5,'k-');

plot(t,y6,'k-');

plot(year,rslt(:,1));

plot(year,rslt(:,3)-0.3);

plot(year,rslt(:,4)-0.6);

plot(year,rslt(:,5)-0.9);

plot(year,rslt(:,6)-1.2);

plot(year,sum(rslt(:,7:8),2)-1.6,'r-');

set(gca,'yTickLabel',[]);

title('EEMD decomposition of GSTA (A_n=0.2; N_e_s_b=100)')

axis([1850 2010 -2.1 0.2]);

xlabel('year');

