% 画出|x|^n+|y|^n=1，按照极坐标有x=r*cos(theta),y=r*sin(theta);
% 因此极坐标方程为r^n*(cos(theta)^n+sin(theta)^n)=1

clear;clc;echo off;

n=5;
dt=1;

figure(1); clf; 
set(gcf,'Color','w'); 
filename = 'unit circle.gif';

% theta=linspace(0,2*pi,400);

x = linspace(-1,1,201);   % 取奇数个点，保证取到0
id = find(~x);
ilist = [(1:10)/10,2:5,10,20,50,100];  % 1:n
for i = ilist
%     r=(1./(abs(cos(theta)).^i+abs(sin(theta)).^i)).^(1/i);
%     polar(theta,r,'r-');
    
    y = (1-abs(x).^i).^(1/i);
    y(id) = 1;
    
%     id = true(1,200);
%     if i<1
%         id = y>3/1000;    %1e-3;
%     end
    
    plot(x,y,'r-','linewidth',2);   % plot(x(id),y(id),'r-','linewidth',2);
    
    hold on
    plot(x,-y,'r-','linewidth',2);
     
    grid on; axis equal;
    axis([-1.5 1.5 -1.2 1.2])
    
%     name=['|x|^' num2str(n) '+|y|^' num2str(n) '=1'];
    
%     text('Interpreter','latex','String','$$|x|^n+|y|^n=1$$',...
%     'Position',[-0.5 -1.9],'FontSize',14,'color','k','fontweight','bold');

    text('Interpreter','latex','String','$$|x|^p+|y|^p=1$$',...
    'Position',[-0.55 0],'FontSize',18,'color','k','fontweight','bold');

    title(['\it\fontname{Times New Roman}\fontsize{20}\color{black}p ' ...
        '\rm= \bf\color{red}',num2str(i)]);
    
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if  i == ilist(1)
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
	end
    
%     print(1,'-dpng',['unit circle ' sprintf('%03d',i) '.png']);
    
    hold off
    
end