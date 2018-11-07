%%%%%%%2018-11-7   Jeff-Zhouxiaolong%%%%%%%%%%%%%
%下面code为计算包络线以及瞬时频率

clear all; clc; close all;
 
fs=400;                                 % 采样频率
N=400;                                  % 数据长度
n=0:1:N-1;
dt=1/fs;
t=n*dt;                                 % 时间序列
A=0.5;                                  % 相位调制幅值
x=(1+0.5*cos(2*pi*5*t)).*cos(2*pi*50*t+A*sin(2*pi*10*t));  % 信号序列
z=hilbert(x');                          % 希尔伯特变换
a=abs(z);                               % 包络线
fnor=instfreq(z);                       % 瞬时频率
fnor=[fnor(1); fnor; fnor(end)];        % 瞬时频率补齐
% 作图
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1), pos(2)-100,pos(3),pos(4)]);
subplot 211; plot(t,x,'k'); hold on;
plot(t,a,'r--','linewidth',2);
title('包络线'); ylabel('幅值'); xlabel(['时间/s' 10 '(a)']);
ylim([-2,2]);
subplot 212; plot(t,fnor*fs,'k'); ylim([43 57]);
title('瞬时频率'); ylabel('频率/Hz');  xlabel(['时间/s' 10 '(b)']);


