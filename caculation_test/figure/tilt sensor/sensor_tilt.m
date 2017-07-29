Ts = [500;500;-2500];
th = 20;
dt = 0.2;
picname1 = ['camera_theta_',num2str(th),'.gif'];
picname2 = ['tilt_pixel_theta_',num2str(th),'.gif'];
theta = th*pi/180;

for fa=0:10:359,
    fai = fa*pi/180;
    
    r3=[sin(theta)*cos(fai);sin(theta)*sin(fai);cos(theta)];
    r1=[
        (sin(theta)*tan(theta)*sin(fai)^2 + cos(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
        -(cos(fai)*sin(fai)*sin(theta)*tan(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
        -(cos(fai)*sin(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
        ];
    r2=cross(r3,r1);
    R=[r1,r2,r3];
    
    cos_alpha = -sin(2*fai)*sin(theta)^2/(sin(2*fai)^2*sin(theta)^4+4*cos(theta)^2)^(1/2);
    alpha = acosd(cos_alpha);
    
    Ks=[-Ts(3),0,Ts(1);0,-Ts(3),Ts(2);0,0,1];
    Tp=[0;0;1]*Ts(3)/r3(3);
    Kp=[-Tp(3),0,Tp(1);0,-Tp(3),Tp(2);0,0,1];
    H=Kp*R'/Ks;
    
    % 2d points: in the local frame
    xs4 = [1000*[0 0 1 1;0 1 1 0];ones(1,4)]+(Ts-r3*Ts(3)/R(3,3))*ones(1,4);
    xp4 = H*xs4;
    xp4 = xp4./(ones(3,1)*xp4(3,:));
    
    dx = xp4(1:2,4);
    dy = xp4(1:2,2);
    lx=norm(dx);
    ly=norm(dy);
    alpha = acosd(dot(dx,dy)/(lx*ly));
    
    
    % draw plane grid
    % vector = (0:0.5:2)*1000;
    % n = length(vector);
    % [xg,yg] = meshgrid(vector);
    % zg = zeros(n);
    % Xg = [xg(:)';yg(:)';zg(:)'];  % in the local frame
    %
    % Xpg = Xg-Tp*ones(1,n*n);   % in the camera frame
    % Xsg = R'*(Xg-Ts*ones(1,n*n));   % in the camera frame
    %
    % zp = reshape(Xpg(3,:),n,n);
    %
    % xs = reshape(Xsg(1,:),n,n);
    % ys = reshape(Xsg(2,:),n,n);
    % zs = reshape(Xsg(3,:),n,n);
    %
    % figure(1);
    % mesh(xg,zp,-yg);
    % hold on;
    % mesh(xs,zs,-ys);
    % axis equal;
    % grid off;
    % set(1,'color','w');
    
    % Draw in the camera frame
    rectangle = [0 0 1 1;0 1 1 0];
    rects = [rectangle;zeros(1,4)]*3000;          % in the local frame
    rects = R'*(rects-Ts*ones(1,4));              % in the camera frame
    rectp = [rectangle-0.3;zeros(1,4)]*3000;      % in the local frame
    rectp = rectp-Tp*ones(1,4);                   % in the camera frame
    
    xsp0 = Ts-r3*Ts(3)/R(3,3); % principal point in the sensor frame
    delta = [-R(2,3);R(1,3);0]*4500;  % direction of intersection line of sensor and image plane
    point1 = R'*(delta-r3*Ts(3)/R(3,3));   % R'*(xsp0+delta-Ts)
    point2 = -R'*(delta+r3*Ts(3)/R(3,3));  % R'*(xsp0-delta-Ts)
    
    figure(1);
    plot3(rects(1,[1:end,1]),rects(3,[1:end,1]),-rects(2,[1:end,1]),'k-','linewidth',2.0);
    hold on;
    plot3(rectp(1,[1:end,1]),rectp(3,[1:end,1]),-rectp(2,[1:end,1]),'k-','linewidth',2.0);
    plot3([point1(1) point2(1)],[point1(3) point2(3)],-[point1(2) point2(2)],'k-','linewidth',1.0);
    axis equal;
    axis off;
    axis vis3d;
    grid off;
    view(45,30);
    
    set(1,'color','w');
    
    % draw frame
    BASE = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]*500;   % in the local frame
    BASE_s = R'*(BASE-Ts*ones(1,6));
    BASE_p = BASE-Tp*ones(1,6);
    
    plot3(BASE(1,:),BASE(3,:),-BASE(2,:),'r-','linewidth',2.0);
    % text(BASE(1,1),BASE(3,1),-BASE(2,1),'\bf\it\fontname{Arial}\fontsize{14}O_c','HorizontalAlignment','center');
    text(BASE(1,2),BASE(3,2),-BASE(2,2),'\bf\it\fontname{Arial}\fontsize{14}X_c','HorizontalAlignment','center');
    text(BASE(1,6),BASE(3,6),-BASE(2,6),'\bf\it\fontname{Arial}\fontsize{14}Z_c','HorizontalAlignment','center');
    text(BASE(1,4),BASE(3,4),-BASE(2,4),'\bf\it\fontname{Arial}\fontsize{14}Y_c','HorizontalAlignment','center');
    
    plot3(BASE_s(1,:),BASE_s(3,:),-BASE_s(2,:),'r-','linewidth',2.0);
    % text(BASE_s(1,1),BASE_s(3,1),-BASE_s(2,1),'\bf\it\fontname{Arial}\fontsize{14}O_s','HorizontalAlignment','center');
    text(BASE_s(1,2),BASE_s(3,2),-BASE_s(2,2),'\bf\it\fontname{Arial}\fontsize{14}X_s','HorizontalAlignment','center');
    text(BASE_s(1,6),BASE_s(3,6),-BASE_s(2,6),'\bf\it\fontname{Arial}\fontsize{14}Z_s','HorizontalAlignment','center');
    text(BASE_s(1,4),BASE_s(3,4),-BASE_s(2,4),'\bf\it\fontname{Arial}\fontsize{14}Y_s','HorizontalAlignment','center');
    
    plot3(BASE_p(1,:),BASE_p(3,:),-BASE_p(2,:),'r-','linewidth',2.0);
    % text(BASE_p(1,1),BASE_p(3,1),-BASE_p(2,1),'\bf\it\fontname{Arial}\fontsize{14}O_p','HorizontalAlignment','center');
    text(BASE_p(1,2),BASE_p(3,2),-BASE_p(2,2),'\bf\it\fontname{Arial}\fontsize{14}X_p','HorizontalAlignment','center');
    text(BASE_p(1,6),BASE_p(3,6),-BASE_p(2,6),'\bf\it\fontname{Arial}\fontsize{14}Z_p','HorizontalAlignment','center');
    text(BASE_p(1,4),BASE_p(3,4),-BASE_p(2,4),'\bf\it\fontname{Arial}\fontsize{14}Y_p','HorizontalAlignment','center');
    
    % draw 3d points: Xs=R*Xc+Ts, Xp=Xc+Tp,so Xc=R'*(Xs-Ts), Xc=Xp-Tp;
    Xs4 = xs4(1:2,:);
    Xs4 = [Xs4;zeros(1,4)];
    Xs4 = R'*(Xs4-Ts*ones(1,4));  % in the camera frame
    
    Xp4 = xp4(1:2,:);
    Xp4 = [Xp4;zeros(1,4)];
    Xp4 = Xp4-Tp*ones(1,4);  % in the camera frame
    
    project_ls = [zeros(3,4);Xs4;zeros(3,4)];
    project_ls = reshape(project_ls,3,12);
    project_lp = [zeros(3,4);Xp4;zeros(3,4)];
    project_lp = reshape(project_lp,3,12);
    
    plot3(project_ls(1,:),project_ls(3,:),-project_ls(2,:),'-','color',0.8*[1 1 1],'linewidth',0.5);
    plot3(project_lp(1,:),project_lp(3,:),-project_lp(2,:),'-','color',0.8*[1 1 1],'linewidth',0.5);
    
    plot3(0,0,0,'k.','markersize',10.0);
    plot3(Xs4(1,:),Xs4(3,:),-Xs4(2,:),'k.','markersize',10.0);
    plot3(Xp4(1,:),Xp4(3,:),-Xp4(2,:),'k.','markersize',10.0);
    plot3(Xs4(1,[1:end,1]),Xs4(3,[1:end,1]),-Xs4(2,[1:end,1]),'b-','linewidth',2.0);
    plot3(Xp4(1,[1:end,1]),Xp4(3,[1:end,1]),-Xp4(2,[1:end,1]),'g-','linewidth',2.0);
    
    text(Xs4(1,2),Xs4(3,2),-Xs4(2,2),'\bf\it\fontname{Arial}\fontsize{14}B','HorizontalAlignment','center');
    text(Xs4(1,3),Xs4(3,3),-Xs4(2,3),'\bf\it\fontname{Arial}\fontsize{14}C','HorizontalAlignment','center');
    text(Xs4(1,4),Xs4(3,4),-Xs4(2,4),'\bf\it\fontname{Arial}\fontsize{14}A','HorizontalAlignment','center');
    
    text(Xp4(1,2),Xp4(3,2),-Xp4(2,2),'\bf\it\fontname{Arial}\fontsize{14}B''','HorizontalAlignment','center');
    text(Xp4(1,3),Xp4(3,3),-Xp4(2,3),'\bf\it\fontname{Arial}\fontsize{14}C''','HorizontalAlignment','center');
    text(Xp4(1,4),Xp4(3,4),-Xp4(2,4),'\bf\it\fontname{Arial}\fontsize{14}A''','HorizontalAlignment','center');
    
    
    set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    axis([-2, 2, -0.1, 5,-2,1]*1000);
    
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    if  fai == 0,
        imwrite(I,map,picname1,'gif','LoopCount',Inf,'DelayTime',dt);
    else
        imwrite(I,map,picname1,'gif','WriteMode','append','DelayTime',dt);
    end;
    hold off;
    
    
    %  A0 = double(subs(A,{theta,fai},[th*pi/180,fa*pi/180]));
    xb = rectangle;
    xg = xp4(1:2,:)/1000;
    
    figure(2);
    plot(xb(1,[1:4 1]),xb(2,[1:4 1]),'b-');
    hold on;
    plot(xb(1,:),xb(2,:),'b.');
    plot(xg(1,[1:4 1]),xg(2,[1:4 1]),'g-');
    plot(xg(1,:),xg(2,:),'g.');
    set(gcf,'color','w');
    set(gca,'ydir','reverse');
    axis equal;
    axis([-0.5 1.5 -0.5 1.5]);
    
    text('Interpreter','latex','String',['$$\alpha = ', num2str(alpha), '$$'],...
        'Position',[0.5 0.5],'FontSize',16,'color','k','fontweight','bold','HorizontalAlignment','center');
    
    title(['$$\theta = ', num2str(th), ';  \varphi = ' num2str(fa), '$$'],'Interpreter','latex',...
        'FontSize',16,'color','k','fontweight','bold');
    
    drawnow;
    frame = getframe(2);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    if  fai == 0,
        imwrite(I,map,picname2,'gif','LoopCount',Inf,'DelayTime',dt);
    else
        imwrite(I,map,picname2,'gif','WriteMode','append','DelayTime',dt);
    end;
    hold off;
end;