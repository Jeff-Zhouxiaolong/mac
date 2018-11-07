%%%%%%%2018-11-7   Jeff-Zhouxiaolong%%%%%%%%%%%%%
%����codeΪ����������Լ�˲ʱƵ��

clear all; clc; close all;
 
fs=400;                                 % ����Ƶ��
N=400;                                  % ���ݳ���
n=0:1:N-1;
dt=1/fs;
t=n*dt;                                 % ʱ������
A=0.5;                                  % ��λ���Ʒ�ֵ
x=(1+0.5*cos(2*pi*5*t)).*cos(2*pi*50*t+A*sin(2*pi*10*t));  % �ź�����
z=hilbert(x');                          % ϣ�����ر任
a=abs(z);                               % ������
fnor=instfreq(z);                       % ˲ʱƵ��
fnor=[fnor(1); fnor; fnor(end)];        % ˲ʱƵ�ʲ���
% ��ͼ
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1), pos(2)-100,pos(3),pos(4)]);
subplot 211; plot(t,x,'k'); hold on;
plot(t,a,'r--','linewidth',2);
title('������'); ylabel('��ֵ'); xlabel(['ʱ��/s' 10 '(a)']);
ylim([-2,2]);
subplot 212; plot(t,fnor*fs,'k'); ylim([43 57]);
title('˲ʱƵ��'); ylabel('Ƶ��/Hz');  xlabel(['ʱ��/s' 10 '(b)']);


