%% dRdom
omega=sym('w%d',[3 1]);
omega=sym(omega,'real');

theta=sqrt(omega(1)^2+omega(2)^2+omega(3)^2);
omega_skew=skew3(omega);
Rc = eye(3) + (sin(theta)/theta)*omega_skew + ((1-cos(theta))/theta^2)*(omega_skew * omega_skew);
dRdom=simplify(jacobian(Rc(:),omega));

syms theta st c1 ct2 w1w1 w2w2 w3w3 w1w2w3 stw1 stw2 stw3 ct2w1w2 ct2w2w3 ct2w1w3 theta2 inv_theta2 cs real
syms Rc11 Rc21 Rc31 Rc12 Rc22 Rc32 Rc13 Rc23 Rc33 real

Rc=subs(Rc,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
Rc=subs(Rc,sin(theta)/theta,st);
Rc=subs(Rc,1-cos(theta),c1);
Rc=subs(Rc,c1/theta^2,ct2);
Rc=subs(Rc,omega(1)^2,w1w1);
Rc=subs(Rc,omega(2)^2,w2w2);
Rc=subs(Rc,omega(3)^2,w3w3);
Rc=subs(Rc,st*omega(1),stw1);
Rc=subs(Rc,st*omega(2),stw2);
Rc=subs(Rc,st*omega(3),stw3);
Rc=subs(Rc,ct2*omega(1)*omega(2),ct2w1w2);
Rc=subs(Rc,ct2*omega(2)*omega(3),ct2w2w3);
Rc=subs(Rc,ct2*omega(1)*omega(3),ct2w1w3);

dRdom=subs(dRdom,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
dRdom=subs(dRdom,sin(theta)/theta,st);
dRdom=subs(dRdom,1-cos(theta),c1);
dRdom=simplify(subs(dRdom,c1/theta^2,ct2));
dRdom=subs(dRdom,omega(1)^2,w1w1);
dRdom=subs(dRdom,omega(2)^2,w2w2);
dRdom=subs(dRdom,omega(3)^2,w3w3);
dRdom=subs(dRdom,omega(1)*omega(2)*omega(3),w1w2w3);

dRdom=subs(dRdom,st*omega(1),stw1);
dRdom=subs(dRdom,st*omega(2),stw2);
dRdom=subs(dRdom,st*omega(3),stw3);
dRdom=subs(dRdom,ct2*omega(1)*omega(2),ct2w1w2);
dRdom=subs(dRdom,ct2*omega(2)*omega(3),ct2w2w3);
dRdom=subs(dRdom,ct2*omega(1)*omega(3),ct2w1w3);
dRdom=subs(dRdom,ct2*omega(1)*omega(3),ct2w1w3);
dRdom=subs(dRdom,Rc(1,1),Rc11);
dRdom=subs(dRdom,Rc(2,1),Rc21);
dRdom=subs(dRdom,Rc(3,1),Rc31);
dRdom=subs(dRdom,Rc(1,2),Rc12);
dRdom=subs(dRdom,Rc(2,2),Rc22);
dRdom=subs(dRdom,Rc(3,2),Rc32);
dRdom=subs(dRdom,Rc(1,3),Rc13);
dRdom=subs(dRdom,Rc(2,3),Rc23);
dRdom=subs(dRdom,Rc(3,3),Rc33);

dRdom=subs(dRdom,theta^2,theta2);
dRdom=subs(dRdom,1/theta^2,inv_theta2);
dRdom=subs(dRdom,cos(theta),cs);

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
% ct2w1=ct2*w1;
% ct2w2=ct2*w2;
% ct2w3=ct2*w3;
% theta2=theta^2;
% inv_theta2=1/theta2;
% cs=cos(theta);

%% domdR

R=sym('r%d%d',[3 3]);
R=sym(R,'real');

% [U,S,V] = svd(R);
%  R = U*V';

xr=(trace(R)-1)/2; 
yrv=[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]/2;
yr=sqrt(yrv(1)^2+yrv(2)^2+yrv(3)^2);
vec=yrv/yr;
theta = atan2(yr,xr);

om = simplify(theta*vec);
domdR=simplify(jacobian(om(:),R(:)));

syms theta xr yr tr2 ant real
domdR = subs(domdR,atan2(((R(1,2)/2 - R(2,1)/2)^2 + (R(1,3)/2 - R(3,1)/2)^2 + (R(2,3)/2 - R(3,2)/2)^2)^(1/2), R(1,1)/2 + R(2,2)/2 + R(3,3)/2 - 1/2),theta);
domdR = subs(domdR,R(1,1)/2 + R(2,2)/2 + R(3,3)/2 - 1/2,xr);
domdR = subs(domdR,((R(1,2)/2 - R(2,1)/2)^2 + (R(1,3)/2 - R(3,1)/2)^2 + (R(2,3)/2 - R(3,2)/2)^2)^(1/2),yr);
domdR = subs(domdR,(R(1,1) + R(2,2) + R(3,3) - 1)^2,tr2);
domdR = subs(domdR,(R(1,2) - R(2,1))^2 + (R(1,3) - R(3,1))^2 + (R(2,3) - R(3,2))^2,ant);
domdR = subs(domdR,(R(1,2)/2 - R(2,1)/2)^2 + (R(1,3)/2 - R(3,1)/2)^2 + (R(2,3)/2 - R(3,2)/2)^2,ant/4);
domdR = simplify(domdR)


% r11=R(1,1);
% r12=R(1,2);
% r13=R(1,3);
% r21=R(2,1);
% r22=R(2,2);
% r23=R(2,3);
% r31=R(3,1);
% r32=R(3,2);
% r33=R(3,3);
% theta = atan2(((r12/2 - r21/2)^2 + (r13/2 - r31/2)^2 + (r23/2 - r32/2)^2)^(1/2), r11/2 + r22/2 + r33/2 - 1/2);
% xr = r11/2 + r22/2 + r33/2 - 1/2;
% yr = ((r12/2 - r21/2)^2 + (r13/2 - r31/2)^2 + (r23/2 - r32/2)^2)^(1/2);
% tr2 = (r11 + r22 + r33 - 1)^2;
% ant = (r12 - r21)^2 + (r13 - r31)^2 + (r23 - r32)^2;



