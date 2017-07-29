% The intrinsic parameters (fc;cc;alpha_c)
fc=sym('fc%d',[2 1]);
fc=sym(fc,'real');
cc=sym('cc%d',[2 1]);
cc=sym(cc,'real');
syms alpha_c real

% The extrinsic parameters (omega;Tc)
omega=sym('w%d',[3 1]);
omega=sym(omega,'real');
Tc=sym('t%d',[3 1]);
Tc=sym(Tc,'real');


% Let X = [X1;X2] be a  model point on the chessboard. 
X=sym('X%d',[2 1]);
X=sym(X,'real');

xm=sym('xm%d',[2 1]);
xm=sym(xm,'real');

% Expression for the rotation matrix based on the Rodrigues formula
theta=sqrt(omega(1)^2+omega(2)^2+omega(3)^2);
omega_skew=skew3(omega);
Rc = eye(3) + (sin(theta)/theta)*omega_skew + ((1-cos(theta))/theta^2)*(omega_skew * omega_skew);

% Let us project now that model point on the image plane according to the intrinsic and extrinsic parameters 
% the intrinsic parameter matrix
KK=[fc(1) alpha_c*fc(1) cc(1); 0 fc(2) cc(2); 0 0 1];
xc = KK*[Rc(:,1),Rc(:,2),Tc]*[X;1];
xp = simplify(xc(1:2)/xc(3));

% calculate the geometric distance in x and y direction
% xm = the pixel positions of an extracted corner
% xp = the pixel positions of the projection of the corresponding model point
dx=xm-xp;

% Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
A=jacobian(dx,[fc;cc;alpha_c]);
B=jacobian(dx,[omega;Tc]);

syms theta st c1 ct2 w1w1 w2w2 w3w3 w1w2w3 ct2w1 ct2w2 ct2w3 theta2 inv_theta2 cs real

A=subs(A,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
B=subs(B,sqrt(omega(1)^2+omega(2)^2+omega(3)^2),theta);
A=subs(A,sin(theta)/theta,st);
B=subs(B,sin(theta)/theta,st);
A=subs(A,1-cos(theta),c1);
B=subs(B,1-cos(theta),c1);
A=subs(A,c1/theta^2,ct2);
B=subs(B,c1/theta^2,ct2);
A=subs(A,omega(1)^2,w1w1);
B=subs(B,omega(1)^2,w1w1);
A=subs(A,omega(2)^2,w2w2);
B=subs(B,omega(2)^2,w2w2);
A=subs(A,omega(3)^2,w3w3);
B=subs(B,omega(3)^2,w3w3);
A=subs(A,omega(1)*omega(2)*omega(3),w1w2w3);
B=subs(B,omega(1)*omega(2)*omega(3),w1w2w3);
A=subs(A,ct2*omega(1),ct2w1);
B=subs(B,ct2*omega(1),ct2w1);
A=subs(A,ct2*omega(2),ct2w2);
B=subs(B,ct2*omega(2),ct2w2);
A=subs(A,ct2*omega(3),ct2w3);
B=subs(B,ct2*omega(3),ct2w3);
A=subs(A,theta^2,theta2);
B=subs(B,theta^2,theta2);
A=subs(A,1/theta^2,inv_theta2);
B=subs(B,1/theta^2,inv_theta2);
A=subs(A,cos(theta),cs);
B=subs(B,cos(theta),cs);

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
% ct2w1=ct2*w1;
% ct2w2=ct2*w2;
% ct2w3=ct2*w3;
% theta2=theta^2;
% inv_theta2=1/theta2;
% cs=cos(theta);