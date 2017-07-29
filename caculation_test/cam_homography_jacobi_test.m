% Test of the Jacobians:
% Homography initialization: no distortion
X = [10*randn(2,1);0];  % chessboard: X=[X1;X2;0]
omega = randn(3,1);
Tc = [10*randn(2,1);100*rand(1)];  % the depth (Zc) of the origin of the chessboard must be positive
fc = 1000*rand(2,1);
cc = 1000*randn(2,1);
kc = zeros(5,1);  % no distortion
alpha_c = 0.01*randn;
[x,dxdom,dxdT,dxdf,dxdc,dxdk,dxdalpha] = project_points2(X,omega,Tc,fc,cc,kc,alpha_c);


%%% Jacobian formular
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
ct2w1=ct2*w1;
ct2w2=ct2*w2;
ct2w3=ct2*w3;
theta2=theta^2;
inv_theta2=1/theta2;
cs=cos(theta);

t1=Tc(1);
t2=Tc(2);
t3=Tc(3);

X1=X(1);
X2=X(2);
% xm1=xm(1);
% xm2=xm(2);


JA(1, 1) = -(t1 + alpha_c*t2 - X2*(alpha_c*(ct2*(w1w1 + w3w3) - 1) - ct2w1*w2 + st*w3) + X1*(alpha_c*(ct2w1*w2 + st*w3) - ct2*(w2w2 + w3w3) + 1))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JA(1, 2) = 0;
JA(1, 3) = -1;
JA(1, 4) = 0;
JA(1, 5) = -(fc1*t2 + X1*fc1*(ct2w1*w2 + st*w3) - X2*fc1*(ct2*(w1w1 + w3w3) - 1))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JA(2, 1) = 0;
JA(2, 2) = -(t2 - X2*(ct2*(w1w1 + w3w3) - 1) + X1*(ct2w1*w2 + st*w3))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JA(2, 3) = 0;
JA(2, 4) = -1;
JA(2, 5) = 0;


JB(1, 1) = ((X2*(st + cs*inv_theta2*w1w1 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 + inv_theta2*st*w1w2w3) + X1*(ct2w3 - 2*ct2w3*inv_theta2*w1w1 - cs*inv_theta2*w1*w2 + inv_theta2*st*w1*w2 + inv_theta2*st*w3*w1w1))*(cc1*t3 + fc1*t1 + X2*(cc1*(ct2w2*w3 + st*w1) + fc1*(ct2w1*w2 - st*w3) - alpha_c*fc1*(ct2*(w1w1 + w3w3) - 1)) + X1*(cc1*(ct2w1*w3 - st*w2) - fc1*(ct2*(w2w2 + w3w3) - 1) + alpha_c*fc1*(ct2w1*w2 + st*w3)) + alpha_c*fc1*t2))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2 - (X1*(cc1*(ct2w3 - 2*ct2w3*inv_theta2*w1w1 - cs*inv_theta2*w1*w2 + inv_theta2*st*w1*w2 + inv_theta2*st*w3*w1w1) + fc1*(2*ct2w1*inv_theta2*(w2w2 + w3w3) - inv_theta2*st*w1*(w2w2 + w3w3)) + alpha_c*fc1*(ct2w2 - 2*ct2w2*inv_theta2*w1w1 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w1w1)) + X2*(cc1*(st + cs*inv_theta2*w1w1 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 + inv_theta2*st*w1w2w3) + fc1*(ct2w2 - 2*ct2w2*inv_theta2*w1w1 - cs*inv_theta2*w1*w3 + inv_theta2*st*w1*w3 + inv_theta2*st*w2*w1w1) - alpha_c*fc1*(2*ct2w1 - 2*ct2w1*inv_theta2*(w1w1 + w3w3) + inv_theta2*st*w1*(w1w1 + w3w3))))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JB(1, 2) = - (X2*(cc1*(ct2w3 - 2*ct2w3*inv_theta2*w2w2 + cs*inv_theta2*w1*w2 - inv_theta2*st*w1*w2 + inv_theta2*st*w3*w2w2) + fc1*(ct2w1 - 2*ct2w1*inv_theta2*w2w2 - cs*inv_theta2*w2*w3 + inv_theta2*st*w2*w3 + inv_theta2*st*w1*w2w2) + alpha_c*fc1*(2*ct2w2*inv_theta2*(w1w1 + w3w3) - inv_theta2*st*w2*(w1w1 + w3w3))) - X1*(cc1*(st + cs*inv_theta2*w2w2 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w2w2 - inv_theta2*st*w1w2w3) + fc1*(2*ct2w2 - 2*ct2w2*inv_theta2*(w2w2 + w3w3) + inv_theta2*st*w2*(w2w2 + w3w3)) - alpha_c*fc1*(ct2w1 - 2*ct2w1*inv_theta2*w2w2 + cs*inv_theta2*w2*w3 - inv_theta2*st*w2*w3 + inv_theta2*st*w1*w2w2)))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1)) - ((X1*(st + cs*inv_theta2*w2w2 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w2w2 - inv_theta2*st*w1w2w3) - X2*(ct2w3 - 2*ct2w3*inv_theta2*w2w2 + cs*inv_theta2*w1*w2 - inv_theta2*st*w1*w2 + inv_theta2*st*w3*w2w2))*(cc1*t3 + fc1*t1 + X2*(cc1*(ct2w2*w3 + st*w1) + fc1*(ct2w1*w2 - st*w3) - alpha_c*fc1*(ct2*(w1w1 + w3w3) - 1)) + X1*(cc1*(ct2w1*w3 - st*w2) - fc1*(ct2*(w2w2 + w3w3) - 1) + alpha_c*fc1*(ct2w1*w2 + st*w3)) + alpha_c*fc1*t2))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2;
JB(1, 3) = (X2*(fc1*(st + cs*inv_theta2*w3w3 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w3w3 - inv_theta2*st*w1w2w3) - cc1*(ct2w2 - 2*ct2w2*inv_theta2*w3w3 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w3w3) + alpha_c*fc1*(2*ct2w3 - 2*ct2w3*inv_theta2*(w1w1 + w3w3) + inv_theta2*st*w3*(w1w1 + w3w3))) - X1*(cc1*(ct2w1 - 2*ct2w1*inv_theta2*w3w3 - cs*inv_theta2*w2*w3 + inv_theta2*st*w2*w3 + inv_theta2*st*w1*w3w3) - fc1*(2*ct2w3 - 2*ct2w3*inv_theta2*(w2w2 + w3w3) + inv_theta2*st*w3*(w2w2 + w3w3)) + alpha_c*fc1*(st + cs*inv_theta2*w3w3 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w3w3 + inv_theta2*st*w1w2w3)))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1)) + ((X1*(ct2w1 - 2*ct2w1*inv_theta2*w3w3 - cs*inv_theta2*w2*w3 + inv_theta2*st*w2*w3 + inv_theta2*st*w1*w3w3) + X2*(ct2w2 - 2*ct2w2*inv_theta2*w3w3 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w3w3))*(cc1*t3 + fc1*t1 + X2*(cc1*(ct2w2*w3 + st*w1) + fc1*(ct2w1*w2 - st*w3) - alpha_c*fc1*(ct2*(w1w1 + w3w3) - 1)) + X1*(cc1*(ct2w1*w3 - st*w2) - fc1*(ct2*(w2w2 + w3w3) - 1) + alpha_c*fc1*(ct2w1*w2 + st*w3)) + alpha_c*fc1*t2))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2;
JB(1, 4) = -fc1/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JB(1, 5) = -(alpha_c*fc1)/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JB(1, 6) = (fc1*(X1 + t1 + X2*alpha_c + alpha_c*t2 + X2*ct2w1*w2 - X1*ct2*w2w2 - X1*ct2*w3w3 - X2*st*w3 + X1*alpha_c*st*w3 + X1*alpha_c*ct2w1*w2 - X2*alpha_c*ct2*w1w1 - X2*alpha_c*ct2*w3w3))/(t3 + X1*ct2w1*w3 + X2*ct2w2*w3 - X1*st*w2 + X2*st*w1)^2;
JB(2, 1) = ((X2*(st + cs*inv_theta2*w1w1 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 + inv_theta2*st*w1w2w3) + X1*(ct2w3 - 2*ct2w3*inv_theta2*w1w1 - cs*inv_theta2*w1*w2 + inv_theta2*st*w1*w2 + inv_theta2*st*w3*w1w1))*(cc2*t3 + fc2*t2 + X1*(cc2*(ct2w1*w3 - st*w2) + fc2*(ct2w1*w2 + st*w3)) - X2*(fc2*(ct2*(w1w1 + w3w3) - 1) - cc2*(ct2w2*w3 + st*w1))))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2 - (X1*(cc2*(ct2w3 - 2*ct2w3*inv_theta2*w1w1 - cs*inv_theta2*w1*w2 + inv_theta2*st*w1*w2 + inv_theta2*st*w3*w1w1) + fc2*(ct2w2 - 2*ct2w2*inv_theta2*w1w1 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w1w1)) + X2*(cc2*(st + cs*inv_theta2*w1w1 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w1w1 + inv_theta2*st*w1w2w3) - fc2*(2*ct2w1 - 2*ct2w1*inv_theta2*(w1w1 + w3w3) + inv_theta2*st*w1*(w1w1 + w3w3))))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JB(2, 2) = - (X2*(cc2*(ct2w3 - 2*ct2w3*inv_theta2*w2w2 + cs*inv_theta2*w1*w2 - inv_theta2*st*w1*w2 + inv_theta2*st*w3*w2w2) + fc2*(2*ct2w2*inv_theta2*(w1w1 + w3w3) - inv_theta2*st*w2*(w1w1 + w3w3))) - X1*(cc2*(st + cs*inv_theta2*w2w2 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w2w2 - inv_theta2*st*w1w2w3) - fc2*(ct2w1 - 2*ct2w1*inv_theta2*w2w2 + cs*inv_theta2*w2*w3 - inv_theta2*st*w2*w3 + inv_theta2*st*w1*w2w2)))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1)) - ((X1*(st + cs*inv_theta2*w2w2 + 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w2w2 - inv_theta2*st*w1w2w3) - X2*(ct2w3 - 2*ct2w3*inv_theta2*w2w2 + cs*inv_theta2*w1*w2 - inv_theta2*st*w1*w2 + inv_theta2*st*w3*w2w2))*(cc2*t3 + fc2*t2 + X1*(cc2*(ct2w1*w3 - st*w2) + fc2*(ct2w1*w2 + st*w3)) - X2*(fc2*(ct2*(w1w1 + w3w3) - 1) - cc2*(ct2w2*w3 + st*w1))))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2;
JB(2, 3) = (X2*(fc2*(2*ct2w3 - 2*ct2w3*inv_theta2*(w1w1 + w3w3) + inv_theta2*st*w3*(w1w1 + w3w3)) - cc2*(ct2w2 - 2*ct2w2*inv_theta2*w3w3 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w3w3)) - X1*(fc2*(st + cs*inv_theta2*w3w3 - 2*ct2*inv_theta2*w1w2w3 - inv_theta2*st*w3w3 + inv_theta2*st*w1w2w3) + cc2*(ct2w1 - 2*ct2w1*inv_theta2*w3w3 - cs*inv_theta2*w2*w3 + inv_theta2*st*w2*w3 + inv_theta2*st*w1*w3w3)))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1)) + ((X1*(ct2w1 - 2*ct2w1*inv_theta2*w3w3 - cs*inv_theta2*w2*w3 + inv_theta2*st*w2*w3 + inv_theta2*st*w1*w3w3) + X2*(ct2w2 - 2*ct2w2*inv_theta2*w3w3 + cs*inv_theta2*w1*w3 - inv_theta2*st*w1*w3 + inv_theta2*st*w2*w3w3))*(cc2*t3 + fc2*t2 + X1*(cc2*(ct2w1*w3 - st*w2) + fc2*(ct2w1*w2 + st*w3)) - X2*(fc2*(ct2*(w1w1 + w3w3) - 1) - cc2*(ct2w2*w3 + st*w1))))/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1))^2;
JB(2, 4) = 0;
JB(2, 5) = -fc2/(t3 + X1*(ct2w1*w3 - st*w2) + X2*(ct2w2*w3 + st*w1));
JB(2, 6) = (fc2*(X2 + t2 + X1*ct2w1*w2 - X2*ct2*w1w1 - X2*ct2*w3w3 + X1*st*w3))/(t3 + X1*ct2w1*w3 + X2*ct2w2*w3 - X1*st*w2 + X2*st*w1)^2;


A=[dxdf,dxdc,dxdalpha]
JA
B=[dxdom,dxdT]
JB
