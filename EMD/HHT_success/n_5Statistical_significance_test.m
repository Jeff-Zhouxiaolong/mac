
clear;

clf;

data=randn(512,1);

rslt=eemd(data,0,1);

imfs=rslt(:,2:8);

[sigline95,logep]=significance(imfs,0.05);

[sigline99,logep]=significance(imfs,0.01);

plot(sigline95(:,1),sigline95(:,2));  %  95 percenta line

hold on

plot(sigline99(:,1),sigline99(:,2),'m-');  % 99 percenta line

plot(logep(:,1),logep(:,2),'r*');  

plot(logep(1,1),logep(1,2),'k*');

grid;

xlabel('LOG2 ( Mean Period )');

ylabel('LOG2 ( Mean Normalized Energy )');

title('Significance test of IMFs of white noise');

axis([0 10 -7 0])

