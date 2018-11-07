clc
clear all
close all
ts = 0.001;
fs = 1/ts;
N = 200;
f = 50;
k = 0:N-1;
t = k*ts;
% �źű任
% ���ۣ�sin�ź�Hilbert�任��Ϊcos�ź�
y = sin(2*pi*f*t);
yh = hilbert(y);    % matlab�����õ��ź��Ǻϳɵĸ��ź�
yi = imag(yh);      % �鲿Ϊ���϶����Hilbert�任
figure
subplot(211)
plot(t, y)
title('ԭʼsin�ź�')
subplot(212)
plot(t, yi)
title('Hilbert�任�ź�')
ylim([-1,1])