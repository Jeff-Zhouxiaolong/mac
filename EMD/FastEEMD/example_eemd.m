% The example shows how to decomposed the signal by EEMD
function [] = example_eemd()   

inputCase = 2;

if(inputCase == 1) % Input Signal
 n = 1024;
 t = 1:n;
 x =  sin((2*pi*t/16)) + sin(2*pi*t/256);
elseif(inputCase == 2) % Blood Flow Velocity
 % example in the paper [22]
 % K. Hu, M. T. Lo, C. K. Peng, Y. Liu, and V. Novak,
 %A Nonlinear Dynamic Approach Reveals a Long-Term Stroke Effect on Cerebral Blood Flow Regulation at Multiple Time Scales
 %, PLoS Computational Biology, 8 (2012), e1002601.
     load BFVL.mat;
     x = BFVL;
end

t = 1:length(x); figure; plot(x,'k','LineWIdth',1); title('input Signal');

fprintf('************ RUN RCADA EEMD ***********\n');

Nstd = 0.4;
NE = 200; % # of ensemble
numImf = 8; % # of imfs
runCEEMD = 0;
maxSift = 10;
typeSpline = 2;
toModifyBC = 2;
randType = 2;
seedNo = 0;
checksignal = 1;

tic;


[imf] = eemd(x,Nstd,NE,numImf); % run EEMD [3]
%[imf] = eemd(x,0,1,numImf); % run EMD [2]
%[imf] = eemd(x,Nstd,NE,numImf,runCEEMD,maxSift,typeSpline,toModifyBC,randType,seedNo,checksignal); % EEMD [3]
%[imf] = eemd(x,Nstd,NE,numImf,1,maxSift,typeSpline,toModifyBC,randType,seedNo,checksignal); % CEEMD [4]
toc;

%imf = imf';

figure;  
nimf = size(imf,1); 
for (m=1:nimf)
  subplot(nimf,1,m); plot(t,imf(m,:),'k','LineWIdth',1.5);
end
subplot(nimf,1,1); title('IMFs');

[1];



