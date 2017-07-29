%%  Xs=Rsc*Xc+Ts, Xo=Xc+To

syms theta fai real
Ts=sym('t%d',[3,1]);
Ts=sym(Ts,'real');
Ks=[-Ts(3),0,Ts(1);0,-Ts(3),Ts(2);0,0,1];
n=[sin(theta)*cos(fai);sin(theta)*sin(fai);cos(theta)];

es1=[1;0;0];
e1=es1+(Ts-es1)*(dot(n,es1)/dot(n,es1-Ts));
e1=simplify(e1./norm(e1));
e2=simplify(cross(n,e1));
R=[e1,e2,n];

To=R'*Ts;
Ko=[-To(3),0,To(1);0,-To(3),To(2);0,0,1];
H=simplify(Ko*R'/Ks);
xs = [0 1 0 1;0 0 1 1;1 1 1 1];
xo=simplify(H*xs);
xo=simplify(xo./(ones(3,1)*xo(3,:)));

K = simplify(Ks*R);

A = [1, 0, -n(1)/n(3);0, 1, -n(2)/n(3);-n(1), -n(2), -n(3)];
uvw = [500;500;2500];

dt = 0.2;
th = 10;
filename = ['tilt_pixel_theta_',num2str(th),'.gif'];
for fa=0:10:359,
    A0 = double(subs(A,{theta,fai},[th*pi/180,fa*pi/180]));
    ts = A0\uvw;
    Hos=double(subs(H,{Ts(1),Ts(2),Ts(3),theta,fai},[ts',th*pi/180,fa*pi/180]));
    xo0=double(subs(xo,{Ts(1),Ts(2),Ts(3),theta,fai},[ts',th*pi/180,fa*pi/180]));
    % xo1=Hos*xs;
    % xo1=xo1./(ones(3,1)*xo1(3,:));
    
    K0=double(subs(K,{Ts(1),Ts(2),Ts(3),theta,fai},[ts',th*pi/180,fa*pi/180]));
    
    dx = norm(xo0(1:2,2));
    dy = norm(xo0(1:2,3));
    alpha = acosd(dot(xo0(1:2,2),xo0(1:2,3))/(dx*dy));
    % alpha = atan2d(xo0(2,3),xo0(1,3));
    
    figure(1);      
    plot(xs(1,[1 2 4 3 1]),xs(2,[1 2 4 3 1]),'g-');
    hold on;
    plot(xs(1,:),xs(2,:),'g.');
    plot(xo0(1,[1 2 4 3 1]),xo0(2,[1 2 4 3 1]),'b-');
    plot(xo0(1,:),xo0(2,:),'b.');
    set(gcf,'color','w');
    set(gca,'ydir','reverse');
    axis equal;
    axis([-0.5 1.5 -0.5 1.5]);
    
    text('Interpreter','latex','String',['$$\alpha = ', num2str(alpha), '$$'],...
    'Position',[0.5 0.5],'FontSize',16,'color','k','fontweight','bold','HorizontalAlignment','center');

    title(['$$\theta = ', num2str(th), ';  \varphi = ' num2str(fa), '$$'],'Interpreter','latex',...
    'FontSize',16,'color','k','fontweight','bold');
    
    drawnow;    
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    if  fa == 0,
		imwrite(I,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
	else
		imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end;
    hold off;
end;