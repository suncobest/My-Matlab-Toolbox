% The intrinsic parameters (fc;cc;alpha_c;kc)
fc=sym('fc%d',[2 1]);
fc=sym(fc,'real');
cc=sym('cc%d',[2 1]);
cc=sym(cc,'real');
syms alpha_c real

kc=sym('kc%d',[5 1]);
kc=sym(kc,'real');

% The extrinsic parameters (omega;Tc)
omega=sym('w%d',[3 1]);
omega=sym(omega,'real');
Tc=sym('t%d',[3 1]);
Tc=sym(Tc,'real');


% Let Xw = [X1;X2;X3] be a  model point in the world reference frame. 
Xw=sym('X%d',[3 1]);
Xw=sym(Xw,'real');

xm=sym('xm%d',[2 1]);
xm=sym(xm,'real');

% Expression for the rotation matrix based on the Rodrigues formula
theta=sqrt(omega(1)^2+omega(2)^2+omega(3)^2);
omega_skew=skew3(omega);
Rc = eye(3) + (sin(theta)/theta)*omega_skew + ((1-cos(theta))/theta^2)*(omega_skew * omega_skew);

% Let Xc be the coordinate vector of P in the camera reference frame. 
% Then Xw and Xc are related to each other through the following rigid motion equation:
Xc = Rc*Xw+Tc;

% Let xn be the normalized (pinhole) image projection: 
xn=Xc(1:2)/Xc(3);

% After including lens distortion, the new normalized point coordinate xd is defined as follows:
r2=xn(1)^2+xn(2)^2;
xd = (1+kc(1)*r2+kc(2)*r2^2+kc(5)*r2^3)*xn + [2*kc(3)*xn(1)*xn(2) + kc(4)*(r2+2*xn(1)^2); kc(3)*(r2+2*xn(2)^2)+2*kc(4)*xn(1)*xn(2)];
% xd = simplify(xd);

% Let us project now that model point on the image plane according to the intrinsic parameters (fc,cc,alpha_c,kc). 
% the intrinsic parameter matrix
KK=[fc(1) alpha_c*fc(1) cc(1); 0 fc(2) cc(2); 0 0 1];
xp=KK*[xd;1];
xp = simplify(xp(1:2));

% calculate the geometric distance in x and y direction
% xm = the pixel positions of an extracted corner
% xp = the pixel positions of the projection of the corresponding model point
dx=xm-xp;

% Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
A=jacobian(dx,[fc;cc;alpha_c;kc]);
B=jacobian(dx,[omega;Tc]);

syms theta st c1 ct2 w1w1 w2w2 w3w3 w1w2w3 stw1 stw2 stw3 ct2w1w2 ct2w2w3 ct2w1w3 theta2 inv_theta2 cs real
syms xd1 xd2 xn1 xn2 Xc1 Xc2 Xc3 Rc11 Rc21 Rc31 Rc12 Rc22 Rc32 Rc13 Rc23 Rc33 lxn2 lxn4 lxn6 real

A=subs(A,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
B=subs(B,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
xd=subs(xd,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
xn=subs(xn,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
Xc=subs(Xc,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
Rc=subs(Rc,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
A=subs(A,sin(theta)/theta,st);
B=subs(B,sin(theta)/theta,st);
xd=subs(xd,sin(theta)/theta,st);
xn=subs(xn,sin(theta)/theta,st);
Xc=subs(Xc,sin(theta)/theta,st);
Rc=subs(Rc,sin(theta)/theta,st);
A=subs(A,1-cos(theta),c1);
B=subs(B,1-cos(theta),c1);
xd=subs(xd,1-cos(theta),c1);
xn=subs(xn,1-cos(theta),c1);
Xc=subs(Xc,1-cos(theta),c1);
Rc=subs(Rc,1-cos(theta),c1);
A=subs(A,c1/theta^2,ct2);
B=subs(B,c1/theta^2,ct2);
xd=subs(xd,c1/theta^2,ct2);
xn=subs(xn,c1/theta^2,ct2);
Xc=subs(Xc,c1/theta^2,ct2);
Rc=subs(Rc,c1/theta^2,ct2);
A=subs(A,omega(1)^2,w1w1);
B=subs(B,omega(1)^2,w1w1);
xd=subs(xd,omega(1)^2,w1w1);
xn=subs(xn,omega(1)^2,w1w1);
Xc=subs(Xc,omega(1)^2,w1w1);
Rc=subs(Rc,omega(1)^2,w1w1);
A=subs(A,omega(2)^2,w2w2);
B=subs(B,omega(2)^2,w2w2);
xd=subs(xd,omega(2)^2,w2w2);
xn=subs(xn,omega(2)^2,w2w2);
Xc=subs(Xc,omega(2)^2,w2w2);
Rc=subs(Rc,omega(2)^2,w2w2);
A=subs(A,omega(3)^2,w3w3);
B=subs(B,omega(3)^2,w3w3);
xd=subs(xd,omega(3)^2,w3w3);
xn=subs(xn,omega(3)^2,w3w3);
Xc=subs(Xc,omega(3)^2,w3w3);
Rc=subs(Rc,omega(3)^2,w3w3);
A=subs(A,omega(1)*omega(2)*omega(3),w1w2w3);
B=subs(B,omega(1)*omega(2)*omega(3),w1w2w3);
xd=subs(xd,omega(1)*omega(2)*omega(3),w1w2w3);
xn=subs(xn,omega(1)*omega(2)*omega(3),w1w2w3);
Xc=subs(Xc,omega(1)*omega(2)*omega(3),w1w2w3);

A=subs(A,st*omega(1),stw1);
B=subs(B,st*omega(1),stw1);
xd=subs(xd,st*omega(1),stw1);
xn=subs(xn,st*omega(1),stw1);
Xc=subs(Xc,st*omega(1),stw1);
Rc=subs(Rc,st*omega(1),stw1);
A=subs(A,st*omega(2),stw2);
B=subs(B,st*omega(2),stw2);
xd=subs(xd,st*omega(2),stw2);
xn=subs(xn,st*omega(2),stw2);
Xc=subs(Xc,st*omega(2),stw2);
Rc=subs(Rc,st*omega(2),stw2);
A=subs(A,st*omega(3),stw3);
B=subs(B,st*omega(3),stw3);
xd=subs(xd,st*omega(3),stw3);
xn=subs(xn,st*omega(3),stw3);
Xc=subs(Xc,st*omega(3),stw3);
Rc=subs(Rc,st*omega(3),stw3);
A=subs(A,ct2*omega(1)*omega(2),ct2w1w2);
B=subs(B,ct2*omega(1)*omega(2),ct2w1w2);
xd=subs(xd,ct2*omega(1)*omega(2),ct2w1w2);
xn=subs(xn,ct2*omega(1)*omega(2),ct2w1w2);
Xc=subs(Xc,ct2*omega(1)*omega(2),ct2w1w2);
Rc=subs(Rc,ct2*omega(1)*omega(2),ct2w1w2);
A=subs(A,ct2*omega(2)*omega(3),ct2w2w3);
B=subs(B,ct2*omega(2)*omega(3),ct2w2w3);
xd=subs(xd,ct2*omega(2)*omega(3),ct2w2w3);
xn=subs(xn,ct2*omega(2)*omega(3),ct2w2w3);
Xc=subs(Xc,ct2*omega(2)*omega(3),ct2w2w3);
Rc=subs(Rc,ct2*omega(2)*omega(3),ct2w2w3);
A=subs(A,ct2*omega(1)*omega(3),ct2w1w3);
B=subs(B,ct2*omega(1)*omega(3),ct2w1w3);
xd=subs(xd,ct2*omega(1)*omega(3),ct2w1w3);
xn=subs(xn,ct2*omega(1)*omega(3),ct2w1w3);
Xc=subs(Xc,ct2*omega(1)*omega(3),ct2w1w3);
Rc=subs(Rc,ct2*omega(1)*omega(3),ct2w1w3);
A=subs(A,theta^2,theta2);
B=subs(B,theta^2,theta2);
A=subs(A,1/theta^2,inv_theta2);
B=subs(B,1/theta^2,inv_theta2);
A=subs(A,cos(theta),cs);
B=subs(B,cos(theta),cs);
A=subs(A,xd(1),xd1);
B=subs(B,xd(1),xd1);
A=subs(A,xd(2),xd2);
B=subs(B,xd(2),xd2);
A=subs(A,xn(1),xn1);
B=subs(B,xn(1),xn1);
A=subs(A,xn(2),xn2);
B=subs(B,xn(2),xn2);
B=subs(B,Xc(1),Xc1);
B=subs(B,Xc(2),Xc2);
B=subs(B,Xc(3),Xc3);
B=subs(B,Rc(1,1),Rc11);
B=subs(B,Rc(2,1),Rc21);
B=subs(B,Rc(3,1),Rc31);
B=subs(B,Rc(1,2),Rc12);
B=subs(B,Rc(2,2),Rc22);
B=subs(B,Rc(3,2),Rc32);
B=subs(B,Rc(1,3),Rc13);
B=subs(B,Rc(2,3),Rc23);
B=subs(B,Rc(3,3),Rc33);

A=subs(A,xn1^2+xn2^2,lxn2);
B=subs(B,xn1^2+xn2^2,lxn2);
A=subs(A,lxn2^2,lxn4);
B=subs(B,lxn2^2,lxn4);
A=subs(A,lxn2^3,lxn6);
B=subs(B,lxn2^3,lxn6);

A=simplify(A)
B=simplify(B)


% w1=omega(1);
% w2=omega(2);
% w3=omega(3);
% theta=sqrt(w1^2+w2^2+w3^2);
% st=sin(theta)/theta;
% c1=1-cos(theta);
% ct2=c1/theta^2;
% w1w1=w1^2;
% w2w2=w2^2;
% w3w3=w3^2;
% w1w2w3=w1*w2*w3;
% stw1=st*w1;
% stw2=st*w2;
% stw3=st*w3;
% ct2w1w2=ct2*w1*w2;
% ct2w2w3=ct2*w2*w3;
% ct2w1w3=ct2*w1*w3;
% theta2=theta^2;
% inv_theta2=1/theta2;
% cs=cos(theta);

% xd1=xd(1);
% xd2=xd(2);
% xn1=xn(1);
% xn2=xn(2);
% lxn2=xn1^2+xn2^2;
% lxn4=lxn2^2;
% lxn6=lxn2^3;
% Xc1=Xc(1);
% Xc2=Xc(2);
% Xc3=Xc(3);
% Rc11=Rc(1,1);
% Rc21=Rc(2,1);
% Rc31=Rc(3,1);
% Rc12=Rc(1,2);
% Rc22=Rc(2,2);
% Rc32=Rc(3,2);
% Rc13=Rc(1,3);
% Rc23=Rc(2,3); 
% Rc33=Rc(3,3);
