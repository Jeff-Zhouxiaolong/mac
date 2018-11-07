
year=gsta(:,1);

inData=gsta(:,2);

rslt=eemd(inData,0,1);

plot(year,rslt(:,2));

hold on;

plot(year,rslt(:,3)-0.3);

plot(year,rslt(:,4)-0.6);

plot(year,rslt(:,5)-0.9);

plot(year,sum(rslt(:,6:8),2)-1.3,'r-');

hold off

set(gca,'yTickLabel',[]);

axis([1850 2010 -1.8 0.3]);

xlabel('year');

%%%%____________________________________________________________________________%%
annotation(figure1,'textbox',...
    [0.822093750000001 0.411304347826087 0.020875 0.0565867066466768],...
    'String',{'R'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.822093750000001 0.492173913043478 0.020875 0.0565867066466768],...
    'String',{'C4'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.822093750000001 0.568695652173913 0.020875 0.0565867066466768],'String',{'C3'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.822093750000001 0.689130434782609 0.020875 0.0565867066466768],...
    'String',{'C2'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.822093750000001 0.863543228385813 0.020875 0.0287356321839081],'String',{'C1'},...
    'FitBoxToText','off',...
    'LineStyle','none');



