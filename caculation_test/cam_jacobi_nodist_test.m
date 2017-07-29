% Test of the Jacobians:
% camera initialization: no distortion
X = 10*randn(3,1); 
omega = randn(3,1);
Tc = [10*randn(2,1);100*rand(1)];  % the depth (Zc) of the origin of the world frame must be positive
fc = 1000*rand(2,1);
cc = 1000*randn(2,1);
kc = zeros(5,1);  % no distortion
alpha_c = 0.01*randn;
[x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points2(X,omega,Tc,fc,cc,kc,alpha_c);


%%% parameters initialization
fc1=fc(1);
fc2=fc(2);
cc1=cc(1);
cc2=cc(2);

w1=omega(1);
w2=omega(2);
w3=omega(3);
theta=sqrt(w1^2+w2^2+w3^2);
st=sin(theta)/theta;
c1=1-cos(theta);
ct2=c1/theta^2;
w1w1=w1^2;
w2w2=w2^2;
w3w3=w3^2;
w1w2w3=w1*w2*w3;
stw1=st*w1;
stw2=st*w2;
stw3=st*w3;
ct2w1w2=ct2*w1*w2;
ct2w2w3=ct2*w2*w3;
ct2w1w3=ct2*w1*w3;
theta2=theta^2;
inv_theta2=1/theta2;
cs=cos(theta);

t1=Tc(1);
t2=Tc(2);
t3=Tc(3);

X1=X(1);
X2=X(2);
X3=X(3);
% xm1=xm(1);
% xm2=xm(2);



xn(1,1) = (t1 - X1*(ct2*(w2w2 + w3w3) - 1) + X2*(ct2w1w2 - stw3) + X3*(ct2w1w3 + stw2))/(t3 - X3*(ct2*(w1w1 + w2w2) - 1) + X1*(ct2w1w3 - stw2) + X2*(ct2w2w3 + stw1));
xn(2,1) = (t2 - X2*(ct2*(w1w1 + w3w3) - 1) + X3*(ct2w2w3 - stw1) + X1*(ct2w1w2 + stw3))/(t3 - X3*(ct2*(w1w1 + w2w2) - 1) + X1*(ct2w1w3 - stw2) + X2*(ct2w2w3 + stw1));

Xc(1,1) = t1 - X1*(ct2*(w2w2 + w3w3) - 1) + X2*(ct2w1w2 - stw3) + X3*(ct2w1w3 + stw2);
Xc(2,1) = t2 - X2*(ct2*(w1w1 + w3w3) - 1) + X3*(ct2w2w3 - stw1) + X1*(ct2w1w2 + stw3);
Xc(3,1) = t3 - X3*(ct2*(w1w1 + w2w2) - 1) + X1*(ct2w1w3 - stw2) + X2*(ct2w2w3 + stw1);

Rc(1, 1) = 1 - ct2*(w2w2 + w3w3);
Rc(1, 2) = ct2w1w2 - stw3;
Rc(1, 3) = ct2w1w3 + stw2;
Rc(2, 1) = ct2w1w2 + stw3;
Rc(2, 2) = 1 - ct2*(w1w1 + w3w3);
Rc(2, 3) = ct2w2w3 - stw1;
Rc(3, 1) = ct2w1w3 - stw2;
Rc(3, 2) = ct2w2w3 + stw1;
Rc(3, 3) = 1 - ct2*(w1w1 + w2w2);

xn1=xn(1);
xn2=xn(2);
lxn2=xn1^2+xn2^2;
lxn4=lxn2^2;
lxn6=lxn2^3;
Xc1=Xc(1);
Xc2=Xc(2);
Xc3=Xc(3);
Rc11=Rc(1,1);
Rc21=Rc(2,1);
Rc31=Rc(3,1);
Rc12=Rc(1,2);
Rc22=Rc(2,2);
Rc32=Rc(3,2);
Rc13=Rc(1,3);
Rc23=Rc(2,3); 
Rc33=Rc(3,3);

%%% Jacobian formular

JA(1, 1) = -(X1*w1w1 + t1*w1w1 + t1*w2w2 + t1*w3w3 + X2*alpha_c*w2w2 + X1*cs*w2w2 + X1*cs*w3w3 + X2*w1*w2 + X3*w1*w3 + alpha_c*t2*w1w1 + alpha_c*t2*w2w2 + alpha_c*t2*w3w3 + X1*alpha_c*w1*w2 + X3*alpha_c*w2*w3 - X2*cs*w1*w2 - X3*cs*w1*w3 - X2*theta*w3*sin(theta) + X3*theta*w2*sin(theta) + X2*alpha_c*cs*w1w1 + X2*alpha_c*cs*w3w3 - X1*alpha_c*cs*w1*w2 - X3*alpha_c*cs*w2*w3 + X1*alpha_c*theta*w3*sin(theta) - X3*alpha_c*theta*w1*sin(theta))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JA(1, 2) = 0;
JA(1, 3) = -1;
JA(1, 4) = 0;
JA(1, 5) = -(fc1*(X2*w2w2 + t2*w1w1 + t2*w2w2 + t2*w3w3 + X2*cs*w1w1 + X2*cs*w3w3 + X1*w1*w2 + X3*w2*w3 - X1*cs*w1*w2 - X3*cs*w2*w3 + X1*theta*w3*sin(theta) - X3*theta*w1*sin(theta)))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JA(2, 1) = 0;
JA(2, 2) = -xn2;
JA(2, 3) = 0;
JA(2, 4) = -1;
JA(2, 5) = 0;


JB(1, 1) = ((X1*w3 + 2*t3*w1 - X3*stw1*w1^2 + X2*theta*sin(theta) - X1*cs*w3 + 2*X3*cs*w1 + X2*cs*w1w1 - X1*stw1*w2 + X2*st*w1w1 + X2*st*w1w2w3 + X1*stw3*w1w1 - X3*stw1*w2w2 - X1*cs*w1*w2)*(X3*cc1*w3w3 + X1*fc1*w1w1 + cc1*t3*w1w1 + cc1*t3*w2w2 + cc1*t3*w3w3 + fc1*t1*w1w1 + fc1*t1*w2w2 + fc1*t1*w3w3 + X1*cc1*w1*w3 + X2*cc1*w2*w3 + X2*fc1*w1*w2 + X3*fc1*w1*w3 + alpha_c*fc1*t2*w1w1 + alpha_c*fc1*t2*w2w2 + alpha_c*fc1*t2*w3w3 + X3*cc1*cs*w1w1 + X3*cc1*cs*w2w2 + X2*alpha_c*fc1*w2w2 + X1*cs*fc1*w2w2 + X1*cs*fc1*w3w3 + X2*alpha_c*cs*fc1*w1w1 + X2*alpha_c*cs*fc1*w3w3 - X1*cc1*cs*w1*w3 - X2*cc1*cs*w2*w3 + X1*alpha_c*fc1*w1*w2 + X3*alpha_c*fc1*w2*w3 - X2*cs*fc1*w1*w2 - X3*cs*fc1*w1*w3 - X1*cc1*theta*w2*sin(theta) + X2*cc1*theta*w1*sin(theta) - X2*fc1*theta*w3*sin(theta) + X3*fc1*theta*w2*sin(theta) + X1*alpha_c*fc1*theta*w3*sin(theta) - X3*alpha_c*fc1*theta*w1*sin(theta) - X1*alpha_c*cs*fc1*w1*w2 - X3*alpha_c*cs*fc1*w2*w3))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta))^2 - (X1*cc1*w3 + 2*X1*fc1*w1 + X2*fc1*w2 + X3*fc1*w3 + 2*cc1*t3*w1 + 2*fc1*t1*w1 - X1*cc1*stw1*w2 + X2*cc1*st*w1w1 + X2*cc1*st*w1w2w3 + X1*cc1*stw3*w1w1 - X3*cc1*stw1*w2w2 - X2*fc1*stw1*w3 + X3*fc1*stw1*w2 + X2*fc1*stw2*w1w1 + X3*fc1*stw3*w1w1 - X1*fc1*stw1*w2w2 - X1*fc1*stw1*w3w3 + 2*alpha_c*fc1*t2*w1 - X3*cc1*stw1*w1^2 + X2*cc1*theta*sin(theta) - X1*cc1*cs*w3 + 2*X3*cc1*cs*w1 + X2*cc1*cs*w1w1 + X1*alpha_c*fc1*w2 - X2*cs*fc1*w2 - X3*cs*fc1*w3 - X1*alpha_c*cs*fc1*w2 + 2*X2*alpha_c*cs*fc1*w1 - X3*alpha_c*cs*fc1*w1w1 + X1*alpha_c*fc1*stw1*w3 - X3*alpha_c*fc1*st*w1w1 + X3*alpha_c*fc1*st*w1w2w3 + X1*alpha_c*fc1*stw2*w1w1 - X2*alpha_c*fc1*stw1*w3w3 - X1*cc1*cs*w1*w2 - X2*cs*fc1*w1*w3 + X3*cs*fc1*w1*w2 - X2*alpha_c*fc1*stw1*w1^2 - X3*alpha_c*fc1*theta*sin(theta) + X1*alpha_c*cs*fc1*w1*w3)/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JB(1, 2) = ((X2*w3 + 2*t3*w2 - X3*stw2*w2^2 - X1*theta*sin(theta) - X2*cs*w3 + 2*X3*cs*w2 - X1*cs*w2w2 + X2*stw1*w2 - X1*st*w2w2 + X1*st*w1w2w3 - X3*stw2*w1w1 + X2*stw3*w2w2 + X2*cs*w1*w2)*(X3*cc1*w3w3 + X1*fc1*w1w1 + cc1*t3*w1w1 + cc1*t3*w2w2 + cc1*t3*w3w3 + fc1*t1*w1w1 + fc1*t1*w2w2 + fc1*t1*w3w3 + X1*cc1*w1*w3 + X2*cc1*w2*w3 + X2*fc1*w1*w2 + X3*fc1*w1*w3 + alpha_c*fc1*t2*w1w1 + alpha_c*fc1*t2*w2w2 + alpha_c*fc1*t2*w3w3 + X3*cc1*cs*w1w1 + X3*cc1*cs*w2w2 + X2*alpha_c*fc1*w2w2 + X1*cs*fc1*w2w2 + X1*cs*fc1*w3w3 + X2*alpha_c*cs*fc1*w1w1 + X2*alpha_c*cs*fc1*w3w3 - X1*cc1*cs*w1*w3 - X2*cc1*cs*w2*w3 + X1*alpha_c*fc1*w1*w2 + X3*alpha_c*fc1*w2*w3 - X2*cs*fc1*w1*w2 - X3*cs*fc1*w1*w3 - X1*cc1*theta*w2*sin(theta) + X2*cc1*theta*w1*sin(theta) - X2*fc1*theta*w3*sin(theta) + X3*fc1*theta*w2*sin(theta) + X1*alpha_c*fc1*theta*w3*sin(theta) - X3*alpha_c*fc1*theta*w1*sin(theta) - X1*alpha_c*cs*fc1*w1*w2 - X3*alpha_c*cs*fc1*w2*w3))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta))^2 - (X2*cc1*w3 + X2*fc1*w1 + 2*cc1*t3*w2 + 2*fc1*t1*w2 + X2*cc1*stw1*w2 - X1*cc1*st*w2w2 + X1*cc1*st*w1w2w3 - X3*cc1*stw2*w1w1 + X2*cc1*stw3*w2w2 - X2*fc1*stw2*w3 + X3*fc1*st*w2w2 + X3*fc1*st*w1w2w3 + X2*fc1*stw1*w2w2 - X1*fc1*stw2*w3w3 + 2*alpha_c*fc1*t2*w2 - X3*cc1*stw2*w2^2 - X1*fc1*stw2*w2^2 - X1*cc1*theta*sin(theta) + X3*fc1*theta*sin(theta) - X2*cc1*cs*w3 + 2*X3*cc1*cs*w2 - X1*cc1*cs*w2w2 + X1*alpha_c*fc1*w1 + 2*X2*alpha_c*fc1*w2 + X3*alpha_c*fc1*w3 + 2*X1*cs*fc1*w2 - X2*cs*fc1*w1 + X3*cs*fc1*w2w2 - X1*alpha_c*cs*fc1*w1 - X3*alpha_c*cs*fc1*w3 + X1*alpha_c*fc1*stw2*w3 - X3*alpha_c*fc1*stw1*w2 - X2*alpha_c*fc1*stw2*w1w1 + X1*alpha_c*fc1*stw1*w2w2 + X3*alpha_c*fc1*stw3*w2w2 - X2*alpha_c*fc1*stw2*w3w3 + X2*cc1*cs*w1*w2 - X2*cs*fc1*w2*w3 + X1*alpha_c*cs*fc1*w2*w3 - X3*alpha_c*cs*fc1*w1*w2)/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JB(1, 3) = ((X1*w1 + X2*w2 + 2*X3*w3 + 2*t3*w3 - X1*cs*w1 - X2*cs*w2 - X1*stw2*w3 + X2*stw1*w3 - X3*stw3*w1w1 - X3*stw3*w2w2 + X1*stw1*w3w3 + X2*stw2*w3w3 - X1*cs*w2*w3 + X2*cs*w1*w3)*(X3*cc1*w3w3 + X1*fc1*w1w1 + cc1*t3*w1w1 + cc1*t3*w2w2 + cc1*t3*w3w3 + fc1*t1*w1w1 + fc1*t1*w2w2 + fc1*t1*w3w3 + X1*cc1*w1*w3 + X2*cc1*w2*w3 + X2*fc1*w1*w2 + X3*fc1*w1*w3 + alpha_c*fc1*t2*w1w1 + alpha_c*fc1*t2*w2w2 + alpha_c*fc1*t2*w3w3 + X3*cc1*cs*w1w1 + X3*cc1*cs*w2w2 + X2*alpha_c*fc1*w2w2 + X1*cs*fc1*w2w2 + X1*cs*fc1*w3w3 + X2*alpha_c*cs*fc1*w1w1 + X2*alpha_c*cs*fc1*w3w3 - X1*cc1*cs*w1*w3 - X2*cc1*cs*w2*w3 + X1*alpha_c*fc1*w1*w2 + X3*alpha_c*fc1*w2*w3 - X2*cs*fc1*w1*w2 - X3*cs*fc1*w1*w3 - X1*cc1*theta*w2*sin(theta) + X2*cc1*theta*w1*sin(theta) - X2*fc1*theta*w3*sin(theta) + X3*fc1*theta*w2*sin(theta) + X1*alpha_c*fc1*theta*w3*sin(theta) - X3*alpha_c*fc1*theta*w1*sin(theta) - X1*alpha_c*cs*fc1*w1*w2 - X3*alpha_c*cs*fc1*w2*w3))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta))^2 - (X1*cc1*w1 + X2*cc1*w2 + 2*X3*cc1*w3 + X3*fc1*w1 + 2*cc1*t3*w3 + 2*fc1*t1*w3 - X1*cc1*stw2*w3 + X2*cc1*stw1*w3 - X3*cc1*stw3*w1w1 - X3*cc1*stw3*w2w2 + X1*cc1*stw1*w3w3 + X2*cc1*stw2*w3w3 + X3*fc1*stw2*w3 - X2*fc1*st*w3w3 + X2*fc1*st*w1w2w3 - X1*fc1*stw3*w2w2 + X3*fc1*stw1*w3w3 + 2*alpha_c*fc1*t2*w3 - X1*fc1*stw3*w3^2 - X2*fc1*theta*sin(theta) - X1*cc1*cs*w1 - X2*cc1*cs*w2 + X3*alpha_c*fc1*w2 + 2*X1*cs*fc1*w3 - X3*cs*fc1*w1 - X2*cs*fc1*w3w3 + 2*X2*alpha_c*cs*fc1*w3 - X3*alpha_c*cs*fc1*w2 + X1*alpha_c*cs*fc1*w3w3 - X3*alpha_c*fc1*stw1*w3 + X1*alpha_c*fc1*st*w3w3 + X1*alpha_c*fc1*st*w1w2w3 - X2*alpha_c*fc1*stw3*w1w1 + X3*alpha_c*fc1*stw2*w3w3 - X1*cc1*cs*w2*w3 + X2*cc1*cs*w1*w3 + X3*cs*fc1*w2*w3 - X2*alpha_c*fc1*stw3*w3^2 + X1*alpha_c*fc1*theta*sin(theta) - X3*alpha_c*cs*fc1*w1*w3)/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JB(1, 4) = -(fc1*(w1w1 + w2w2 + w3w3))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JB(1, 5) = -(alpha_c*fc1*(w1w1 + w2w2 + w3w3))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta));
JB(1, 6) = (fc1*(w1w1 + w2w2 + w3w3)*(X1*w1w1 + t1*w1w1 + t1*w2w2 + t1*w3w3 + X2*alpha_c*w2w2 + X1*cs*w2w2 + X1*cs*w3w3 + X2*w1*w2 + X3*w1*w3 + alpha_c*t2*w1w1 + alpha_c*t2*w2w2 + alpha_c*t2*w3w3 + X1*alpha_c*w1*w2 + X3*alpha_c*w2*w3 - X2*cs*w1*w2 - X3*cs*w1*w3 - X2*theta*w3*sin(theta) + X3*theta*w2*sin(theta) + X2*alpha_c*cs*w1w1 + X2*alpha_c*cs*w3w3 - X1*alpha_c*cs*w1*w2 - X3*alpha_c*cs*w2*w3 + X1*alpha_c*theta*w3*sin(theta) - X3*alpha_c*theta*w1*sin(theta)))/(X3*w3w3 + t3*w1w1 + t3*w2w2 + t3*w3w3 + X3*cs*w1w1 + X3*cs*w2w2 + X1*w1*w3 + X2*w2*w3 - X1*cs*w1*w3 - X2*cs*w2*w3 - X1*theta*w2*sin(theta) + X2*theta*w1*sin(theta))^2;
JB(2, 1) = (fc2*(X3*(st + cs*inv_theta2*w1w1 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 - inv_theta2*st*w1w2w3) + X2*(2*ct2*w1 + inv_theta2*stw1*(w1w1 + w3w3) - 2*ct2*inv_theta2*w1*(w1w1 + w3w3)) - X1*(ct2*w2 - inv_theta2*stw1*w3 + inv_theta2*stw2*w1w1 + cs*inv_theta2*w1*w3 - 2*ct2*inv_theta2*w2*w1w1)))/Xc3 + (fc2*xn2*(X2*(st + cs*inv_theta2*w1w1 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 + inv_theta2*st*w1w2w3) - X3*(2*ct2*w1 + inv_theta2*stw1*(w1w1 + w2w2) - 2*ct2*inv_theta2*w1*(w1w1 + w2w2)) + X1*(ct2*w3 + inv_theta2*stw1*w2 + inv_theta2*stw3*w1w1 - cs*inv_theta2*w1*w2 - 2*ct2*inv_theta2*w3*w1w1)))/Xc3;
JB(2, 2) = - (fc2*(X1*(ct2*w1 - inv_theta2*stw2*w3 + inv_theta2*stw1*w2w2 + cs*inv_theta2*w2*w3 - 2*ct2*inv_theta2*w1*w2w2) - X2*(inv_theta2*stw2*(w1w1 + w3w3) - 2*ct2*inv_theta2*w2*(w1w1 + w3w3)) + X3*(ct2*w3 + inv_theta2*stw1*w2 + inv_theta2*stw3*w2w2 - cs*inv_theta2*w1*w2 - 2*ct2*inv_theta2*w3*w2w2)))/Xc3 - (fc2*xn2*(X1*(st + cs*inv_theta2*w2w2 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w2w2 - inv_theta2*st*w1w2w3) + X3*(2*ct2*w2 + inv_theta2*stw2*(w1w1 + w2w2) - 2*ct2*inv_theta2*w2*(w1w1 + w2w2)) - X2*(ct2*w3 - inv_theta2*stw1*w2 + inv_theta2*stw3*w2w2 + cs*inv_theta2*w1*w2 - 2*ct2*inv_theta2*w3*w2w2)))/Xc3;
JB(2, 3) = (fc2*xn2*(X1*(ct2*w1 + inv_theta2*stw2*w3 + inv_theta2*stw1*w3w3 - cs*inv_theta2*w2*w3 - 2*ct2*inv_theta2*w1*w3w3) - X3*(inv_theta2*stw3*(w1w1 + w2w2) - 2*ct2*inv_theta2*w3*(w1w1 + w2w2)) + X2*(ct2*w2 - inv_theta2*stw1*w3 + inv_theta2*stw2*w3w3 + cs*inv_theta2*w1*w3 - 2*ct2*inv_theta2*w2*w3w3)))/Xc3 - (fc2*(X1*(st + cs*inv_theta2*w3w3 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w3w3 + inv_theta2*st*w1w2w3) - X2*(2*ct2*w3 + inv_theta2*stw3*(w1w1 + w3w3) - 2*ct2*inv_theta2*w3*(w1w1 + w3w3)) + X3*(ct2*w2 + inv_theta2*stw1*w3 + inv_theta2*stw2*w3w3 - cs*inv_theta2*w1*w3 - 2*ct2*inv_theta2*w2*w3w3)))/Xc3;
JB(2, 4) = 0;
JB(2, 5) = -fc2/Xc3;
JB(2, 6) = (fc2*xn2)/Xc3;



A=[dxdf,dxdc,dxdalpha]
JA
B=[dxdom,dxdT]
JB

