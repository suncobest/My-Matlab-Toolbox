%%  x=diag(R,1)*E_43(-1/f)*diag(R^T,1)*X

% syms theta fai real
% n=[sin(theta)*cos(fai);sin(theta)*sin(fai);cos(theta)];

filename = 'Imaging process.gif';

e1 = -60;
f1 = 85;
h1 = e1/2;

% e2 = 80;
% f2 = 105;
% 
% d = -100;
% T = [0;0;d];

cube = [0 1 1 0 0 1 1 0
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1]-ones(3,8)/2;
abcsize = 50*[3;2;1];
om = 0.1*[-2;-1;3 ];  %  [0;0;0];% rand(3,1);
Rw = rodrigues(om);

%% depth of object
So = 250;

Tw = [0;0;So];        % rand(3,1);
Xw = Rw*diag(abcsize)*cube+Tw*ones(1,8);
Xwl = [Xw(:,1:4),Xw(:,1),Xw(:,5:8),Xw(:,5:6),Xw(:,2:3),Xw(:,7:8),Xw(:,4)]; % cube edge

% principle plane
BASE = abs(f1) *([0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]); 
rectP1 = abs(f1)*[1 -1 -1 1;1 1 -1 -1;0 0 0 0]/2 + h1 * repmat([0;0;1],1,4);
rectP2 = abs(f1)*[1 -1 -1 1;1 1 -1 -1;0 0 0 0]/2 + (h1 - e1) * repmat([0;0;1],1,4);
rectIm = max(abcsize)*[1 -1 -1 1;1 1 -1 -1;0 0 0 0] +(h1- (e1+3*f1)) * repmat([0;0;1],1,4);

%% rotate the lens
da = 6;
theta = 20 *pi/180;
% fai = (0:da:359) *pi/180;
dt = 0.2;
for fai = (0:da:359) *pi/180,    
    R = [cos(fai),-sin(fai),0; sin(fai),cos(fai),0; 0,0,1]*[cos(theta),0,sin(theta); 0,1,0; -sin(theta),0,cos(theta)];
    H1 = [1,0,0,0;0,1,0,0;0,0,1+(e1-h1)/f1,(h1-e1)*h1/f1-e1;0,0,-1/f1,1+h1/f1];  % H=E_34(-e)*E_43(-1/f)，透镜组正则透视矩阵  
    Hsing = [R,zeros(3,1); 0,0,0,1]*H1*[R',zeros(3,1);0,0,0,1]; 
    
    % H2 = [1,0,0,0;0,1,0,0;0,0,1+e2/f2,-e2;0,0,-1/f2,1];
    % Hdou = [R,T; 0,0,0,1]*H2*[R',-R'*T;0,0,0,1]*H1;
    
    BASE_i = R*BASE;
    rect1i = R*rectP1;
    rect2i = R*rectP2;
    
    % image
    Xi=Hsing*[Xw;ones(1,8)];
    Xi=Xi./(ones(4,1)*Xi(4,:));
    Xi(4,:)=[];
    Xil = [Xi(:,1:4),Xi(:,1),Xi(:,5:8),Xi(:,5:6),Xi(:,2:3),Xi(:,7:8),Xi(:,4)];  % edge
    
    % plane coordinate(projective geometry) P'*X=0, P1*X1=0; if X1=H*X, we have P1=inv(H')*P1
    % plane1 = [0;0;1;-h1];
    % plane2 = [0;0;1;e1-h1];
    % plane1i = R*P_plane1; % inv(R')=R
    % plane2i = R*P_plane2;
    
    plane1 = [R(:,3);-h1];     % 物方主平面
    plane2 = [R(:,3);e1-h1];   % 像方主平面
    % 物点到像方主平面的垂线与像方主平面的交点
    lamda2 = -plane2'* [Xw;ones(1,8)];
    Xw_p2 = Xw+R(:,3)*lamda2;
    % 像点到物方主平面的垂线与物方主平面的交点
    lamda1 = -plane1'* [Xi;ones(1,8)];
    Xi_p1 = Xi+R(:,3)*lamda1;
    
    
    % principle point and focus
    Oo = R*[0;0;h1];
    Oi = R*[0;0;h1-e1];  
    Fo = (f1+h1)*R(:,3);
    Fi = (h1-(e1+f1))*R(:,3);
    
    
    figure(1);
    hold off;
    plot3(BASE_i(1,:),BASE_i(3,:),-BASE_i(2,:),'r-','linewidth',2.0);
    hold on;
    axis equal;
    set(1,'color','w');
    text(BASE_i(1,2),BASE_i(3,2),-BASE_i(2,2),'\bf\it\fontname{Arial}\fontsize{14}x','HorizontalAlignment','center');
    text(BASE_i(1,6),BASE_i(3,6),-BASE_i(2,6),'\bf\it\fontname{Arial}\fontsize{14}z','HorizontalAlignment','center');
    text(BASE_i(1,4),BASE_i(3,4),-BASE_i(2,4),'\bf\it\fontname{Arial}\fontsize{14}y','HorizontalAlignment','center');
    
    plot3(rect1i(1,[1:4,1]),rect1i(3,[1:4,1]),-rect1i(2,[1:4,1]),'b-','linewidth',2.0);
    plot3(rect2i(1,[1:4,1]),rect2i(3,[1:4,1]),-rect2i(2,[1:4,1]),'b-','linewidth',2.0);
    text(rect1i(1,1),rect1i(3,1),-rect1i(2,1),'\bf\it\fontname{Arial}\fontsize{14}P_o','HorizontalAlignment','center');
    text(rect2i(1,1),rect2i(3,1),-rect2i(2,1),'\bf\it\fontname{Arial}\fontsize{14}P_i','HorizontalAlignment','center');
    plot3(Oo(1),Oo(3),-Oo(2),'b+','markersize',8.0,'linewidth',1.5);
    plot3(Oi(1),Oi(3),-Oi(2),'b+','markersize',8.0,'linewidth',1.5);
    plot3(Fo(1),Fo(3),-Fo(2),'b+','markersize',8.0,'linewidth',1.5);
    plot3(Fi(1),Fi(3),-Fi(2),'b+','markersize',8.0,'linewidth',1.5);
    plot3([Oi(1);Oo(1)],[Oi(3);Oo(3)],-[Oi(2);Oo(2)],'-','color',0.8*[1 1 1],'linewidth',1.0);    % optical axis
    
    
    % scene
    plot3(Xw(1,:),Xw(3,:),-Xw(2,:),'r.','markersize',10.0);
    plot3(Xwl(1,:),Xwl(3,:),-Xwl(2,:),'g-','linewidth',2.0);
    
    plot3(Xi(1,:),Xi(3,:),-Xi(2,:),'r.','markersize',10.0);
    plot3(Xil(1,:),Xil(3,:),-Xil(2,:),'g-','linewidth',2.0);
    
    plot3([Xw(1,:);Oo(1)*ones(1,8)],[Xw(3,:);Oo(3)*ones(1,8)],-[Xw(2,:);Oo(2)*ones(1,8)],'-','color',0.8*[1 1 1],'linewidth',0.5);
    plot3([Xi(1,:);Oi(1)*ones(1,8)],[Xi(3,:);Oi(3)*ones(1,8)],-[Xi(2,:);Oi(2)*ones(1,8)],'-','color',0.8*[1 1 1],'linewidth',0.5);
    
    plot3(Xw_p2(1,:),Xw_p2(3,:),-Xw_p2(2,:),'r.','markersize',10.0);
    plot3(Xi_p1(1,:),Xi_p1(3,:),-Xi_p1(2,:),'r.','markersize',10.0);
    
    XX = [Xw(1,:);Xw_p2(1,:);Xi(1,:)];
    YY = [Xw(3,:);Xw_p2(3,:);Xi(3,:)];
    ZZ = -[Xw(2,:);Xw_p2(2,:);Xi(2,:)];
    plot3(XX,YY,ZZ,'-','color',0.8*[1 1 1],'linewidth',0.5);
    
    XX = [Xw(1,:);Xi_p1(1,:);Xi(1,:)];
    YY = [Xw(3,:);Xi_p1(3,:);Xi(3,:)];
    ZZ = -[Xw(2,:);Xi_p1(2,:);Xi(2,:)];
    plot3(XX,YY,ZZ,'-','color',0.8*[1 1 1],'linewidth',0.5);
    
    % draw 2d image
    % image plane: Z=h1-(e1+3*f1);
    lamda_xi = ((h1-(e1+3*f1))*ones(1,8)-Xi(3,:))./(Xi(3,:)-Oi(3)*ones(1,8));
    xi2d = Xi+(Xi-Oi*ones(1,8)).*lamda_xi(ones(3,1),:);
    xi2dl = [xi2d(:,1:4),xi2d(:,1),xi2d(:,5:8),xi2d(:,5:6),xi2d(:,2:3),xi2d(:,7:8),xi2d(:,4)];
    
    plot3([xi2d(1,:);Oi(1)*ones(1,8)],[xi2d(3,:);Oi(3)*ones(1,8)],-[xi2d(2,:);Oi(2)*ones(1,8)],'-','color',0.8*[1 1 1],'linewidth',0.5);
    plot3(rectIm(1,[1:4,1]),rectIm(3,[1:4,1]),-rectIm(2,[1:4,1]),'b-','linewidth',2.0);
    plot3(xi2d(1,:),xi2d(3,:),-xi2d(2,:),'k.','markersize',8.0);
    plot3(xi2dl(1,:),xi2dl(3,:),-xi2dl(2,:),'-','color',[1 0 1],'linewidth',1.5);
    
    
    set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    axis([-4, 4, -5.5, 4,-2.5,3]*So/5);  
    view(120,20);
    axis vis3d;
    grid off;
    axis off;
%     axis tight;
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    
    if fai == 0,
		imwrite(I,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
	else
		imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end;  
    
end;


