syms xs ys theta fai real
kesi=[xs;ys;0];

R=sym('r%d%d',[3,3]);
R=sym(R,'real');
r3=R(:,3);

Ts=sym('t%d',[3,1]);
Ts=sym(Ts,'real');
Ks=[-Ts(3),0,Ts(1);0,-Ts(3),Ts(2);0,0,1];
Tp=[0;0;1]*Ts(3)/r3(3);
Kp=[-Tp(3),0,Tp(1);0,-Tp(3),Tp(2);0,0,1];

H=simplify(Kp*R'/Ks)

yita = -(Ts(3)/(R(3,3)*(r3'*(Ts-kesi))))*R'*(Ts-kesi)+Tp;
yita=simplify(yita)

invHt=Kp'\R'*Ks'   % invHt=inv(H')
lm=invHt*[0 1;1 0;Ts(3)*R(2,3)/R(3,3)-Ts(2) Ts(3)*R(1,3)/R(3,3)-Ts(1)];  % 两条直线的单应
lm=simplify(lm)

r3=[sin(theta)*cos(fai);sin(theta)*sin(fai);cos(theta)];
r1=[
 (sin(theta)*tan(theta)*sin(fai)^2 + cos(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
      -(cos(fai)*sin(fai)*sin(theta)*tan(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
                          -(cos(fai)*sin(theta))/(sin(fai)^2*tan(theta)^2 + 1)^(1/2)
];
r2=simplify(cross(r3,r1));
R=[r1,r2,r3];
l=[R(2,1)*R(3,3)-R(2,3)*R(3,1);R(2,2)*R(3,3)-R(2,3)*R(3,2);0];
l=simplify(l)
m=[R(1,1)*R(3,3)-R(1,3)*R(3,1);R(1,2)*R(3,3)-R(1,3)*R(3,2);0];
m=simplify(m)
cos_alpha = dot(l,m)/(norm(l)*norm(m))

%%

Tp=[0;0;1]*Ts(3)/r3(3);
Kp=[-Tp(3),0,Tp(1);0,-Tp(3),Tp(2);0,0,1];
H=Kp*R'/Ks;
xs4 = [0 0 1 1;0 1 1 0;1 1 1 1]+(Ts-r3*Ts(3)/R(3,3))*ones(1,4);
xp4 = H*xs4;
xp4 = simplify(xp4./(ones(3,1)*xp4(3,:)))

%%
syms t1 t2 t3 theta fai real
t = [t1;t2;t3];      % the origin of the camera center in the sensor coordinate system Xs0 
k = [-t3 0 t1;0 -t3 t2;0 0 1];
R = [cos(theta)+(1-cos(theta))*cos(fai)^2,(1-cos(theta))*cos(fai)*sin(fai),sin(theta)*sin(fai);(1-cos(theta))*cos(fai)*sin(fai),cos(theta)+(1-cos(theta))*sin(fai)^2,-sin(theta)*cos(fai);-sin(theta)*sin(fai),sin(theta)*cos(fai),cos(theta)];
T = simplify(R*t);   % the origin of the camera center in the sensor coordinate system Xs1 
K = [-T(3) 0 T(1);0 -T(3) T(2);0 0 1];

H = simplify(K*R/k);          
xs0 = [0 1 0 1;0 0 1 1;1 1 1 1];
xs = H*xs0;

% Xs0=Xc+t, Xs1=R*Xc+T=R(Xc+t)=R*Xs0, T=R*t
% xs0=k*xn, xs1=K*R*xn, xs1=K*R*inv(k)*xs0
% H=K*R*inv(k), xs1=H*xs0


figure(1);
hold on;
set(gcf,'color','w');
set(gca,'ydir','reverse');
plot(xs0(1,[1 2 4 3 1]),xs0(2,[1 2 4 3 1]),'g-');
plot(xs0(1,:),xs0(2,:),'g+');
plot(xs1(1,[1 2 4 3 1]),xs1(2,[1 2 4 3 1]),'b-');
plot(xs1(1,:),xs1(2,:),'r+');
axis([-1 2 -1 2]);
axis equal;

hold off;


%%
syms t1 t2 t3 theta fai real
n=[sin(theta)*cos(fai);sin(theta)*sin(fai);cos(theta)];

a1=(t2*sin(theta)*sin(fai)+t3*cos(theta))/(dot([t1-1;t2;t3],n));
a2=-t2*sin(theta)*cos(fai)/(dot([t1-1;t2;t3],n));
a3=-t3*sin(theta)*cos(fai)/(dot([t1-1;t2;t3],n));

b1=-t1*sin(theta)*sin(fai)/(dot([t1;t2-1;t3],n));
b2=(t1*sin(theta)*cos(fai)+t3*cos(theta))/(dot([t1;t2-1;t3],n));
b3=-t3*sin(theta)*sin(fai)/(dot([t1;t2-1;t3],n));

c1=((t2-t1)*sin(theta)*sin(fai)+t3*cos(theta))/(dot([t1-1;t2-1;t3],n));
c2=((t1-t2)*sin(theta)*cos(fai)+t3*cos(theta))/(dot([t1-1;t2-1;t3],n));
c3=-t3*sin(theta)*(sin(fai)+cos(fai))/(dot([t1-1;t2-1;t3],n));

e1=[a1;a2;a3];
e1=e1*dot([t1-1;t2;t3],n);
e2=cross(n,e1);
e2=simplify(e2/norm(e1));
e1=simplify(e1/norm(e1));
R=[e1,e2,n];
u=[0 b1 c1 a1;0 b2 c2 a2;0 b3 c3 a3];
v=simplify(R'*u);
A=v(:,4);
B=v(:,2);
C=v(:,3);
AA=[((t2*sin(theta)+t3*cos(theta)*sin(fai))^2+t3^2*cos(fai)^2);0;0]/(dot([t1-1;t2;t3],n)*sqrt((t2*sin(theta)*sin(fai)+t3*cos(theta))^2+(t2^2+t3^2)*sin(theta)^2*cos(fai)^2));
BB=[-sin(theta)*(t1*t2*sin(theta)+t2*t3*cos(theta)*cos(fai)+t1*t3*cos(theta)*sin(fai)-t3^2*sin(theta)*sin(fai)*cos(fai));t3*(t1*sin(theta)*cos(fai)+t2*sin(theta)*sin(fai)+t3*cos(theta));0]/(dot([t1;t2-1;t3],n)*sqrt((t2*sin(theta)*sin(fai)+t3*cos(theta))^2+(t2^2+t3^2)*sin(theta)^2*cos(fai)^2));
CC=[(t2*sin(theta)+t3*cos(theta)*sin(fai))^2+t3^2*cos(fai)^2-sin(theta)*(t1*t2*sin(theta)+t2*t3*cos(theta)*cos(fai)+t1*t3*cos(theta)*sin(fai)-t3^2*sin(theta)*sin(fai)*cos(fai));t3*(t1*sin(theta)*cos(fai)+t2*sin(theta)*sin(fai)+t3*cos(theta));0]/(dot([t1-1;t2-1;t3],n)*sqrt((t2*sin(theta)*sin(fai)+t3*cos(theta))^2+(t2^2+t3^2)*sin(theta)^2*cos(fai)^2));
simplify(A-AA)
simplify(B-BB)
simplify(C-CC)

%%
syms xs ys t1 t2 t3 theta fai real
R=sym('r%d%d',3);
R=sym(R,'real');
y=(dot(R(:,3),[xs;ys;0])/dot(R(:,3),([xs;ys;0]-[t1;t2;t3])))*([t1;t2;t3]-[xs;ys;0])+[xs;ys;0];
x=simplify(R'*y);

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
th = 10;
fa = 30;

A0 = double(subs(A,{theta,fai},[th*pi/180,fa*pi/180]));
ts = A0\uvw;
Hos=double(subs(H,{Ts(1),Ts(2),Ts(3),theta,fai},[ts',th*pi/180,fa*pi/180]));
xo0=double(subs(xo,{Ts(1),Ts(2),Ts(3),theta,fai},[ts',th*pi/180,fa*pi/180]));


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
