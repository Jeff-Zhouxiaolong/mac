clc
clear all
close all
ts = 0.001;
fs = 1/ts;
N = 200;
f = 50;
k = 0:N-1;
t = k*ts;
% 信号变换
% 结论：sin信号Hilbert变换后为cos信号
y = sin(2*pi*f*t);
yh = hilbert(y);    % matlab函数得到信号是合成的复信号
yi = imag(yh);      % 虚部为书上定义的Hilbert变换
figure
subplot(211)
plot(t, y)
title('原始sin信号')
subplot(212)
plot(t, yi)
title('Hilbert变换信号')
ylim([-1,1])