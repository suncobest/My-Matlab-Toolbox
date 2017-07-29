% \copyright: zjliu
% Author's email: zjliu2001@163.com 

clc;close all;clear;
axis([-4,4,-4,4]);hold on;
title('摆线的绘制');

set(gcf,'DoubleBuffer','on');
axis square;
tq=linspace(0,pi*2,200);
plot(2*exp(i*tq),'k');
z=-2;
xx=z+exp(i*tq)/4;
hc=plot(xx,'r');
hp=plot(real(xx(1)),imag(xx(1)),'b*');
ht=plot(real(xx(1)),imag(xx(1)),'b');   % 摆线
t=0;dt=0.02;
zk=[xx(1)];
omega=20;  % 转速
while t<8;
    t=t+dt;
    dp=t*omega;
    z=2*exp(i*[pi*(1-t)]);
    xx=z+exp(i*(tq+dp))/4;
    zk=[zk,xx(1)];
    set(hc,'XData',real(xx),'YData',imag(xx));
    set(hp,'XData',real(xx(1)),'YData',imag(xx(1)));
    set(ht,'XData',real(zk),'YData',imag(zk));
    pause(0.1);
end