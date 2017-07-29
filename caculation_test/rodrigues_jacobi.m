%% dRdom
omega = randn(3,1);

dom = randn(3,1)/1000000;
[R1,dR1] = rodrigues(omega);  % dR1 is the derivative of R1 wrt om(矩阵对向量求导)
R2 = rodrigues(omega+dom); 
R2a = R1 + reshape(dR1 * dom,3,3);  
gain = norm(R2 - R1)/norm(R2 - R2a) 


w1=omega(1);
w2=omega(2);
w3=omega(3);
theta=sqrt(w1^2+w2^2+w3^2);

omega_skew=skew3(omega);
Rc = eye(3) + (sin(theta)/theta)*omega_skew + ((1-cos(theta))/theta^2)*(omega_skew * omega_skew);

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


dRdom(1, 1) = inv_theta2*w1*(w2w2 + w3w3)*(2*ct2 - st);
dRdom(1, 2) = 2*ct2*inv_theta2*w2*(w2w2 + w3w3) - inv_theta2*stw2*(w2w2 + w3w3) - 2*ct2*w2;
dRdom(1, 3) = 2*ct2*inv_theta2*w3*(w2w2 + w3w3) - inv_theta2*stw3*(w2w2 + w3w3) - 2*ct2*w3;
dRdom(2, 1) = inv_theta2*(stw2*w1w1 - stw1*w3 + ct2*theta2*w2 + cs*w1*w3 - 2*ct2*w2*w1w1);
dRdom(2, 2) = inv_theta2*(stw1*w2w2 - stw2*w3 + ct2*theta2*w1 + cs*w2*w3 - 2*ct2*w1*w2w2);
dRdom(2, 3) = inv_theta2*(cs*w3w3 - 2*ct2*w1w2w3 + st*theta2 - st*w3w3 + st*w1w2w3);
dRdom(3, 1) = inv_theta2*(stw1*w2 + stw3*w1w1 + ct2*theta2*w3 - cs*w1*w2 - 2*ct2*w3*w1w1);
dRdom(3, 2) = -inv_theta2*(cs*w2w2 + 2*ct2*w1w2w3 + st*theta2 - st*w2w2 - st*w1w2w3);
dRdom(3, 3) = inv_theta2*(stw2*w3 + stw1*w3w3 + ct2*theta2*w1 - cs*w2*w3 - 2*ct2*w1*w3w3);
dRdom(4, 1) = inv_theta2*(stw1*w3 + stw2*w1w1 + ct2*theta2*w2 - cs*w1*w3 - 2*ct2*w2*w1w1);
dRdom(4, 2) = inv_theta2*(stw2*w3 + stw1*w2w2 + ct2*theta2*w1 - cs*w2*w3 - 2*ct2*w1*w2w2);
dRdom(4, 3) = -inv_theta2*(cs*w3w3 + 2*ct2*w1w2w3 + st*theta2 - st*w3w3 - st*w1w2w3);
dRdom(5, 1) = 2*ct2*inv_theta2*w1*(w1w1 + w3w3) - inv_theta2*stw1*(w1w1 + w3w3) - 2*ct2*w1;
dRdom(5, 2) = inv_theta2*w2*(w1w1 + w3w3)*(2*ct2 - st);
dRdom(5, 3) = 2*ct2*inv_theta2*w3*(w1w1 + w3w3) - inv_theta2*stw3*(w1w1 + w3w3) - 2*ct2*w3;
dRdom(6, 1) = inv_theta2*(cs*w1w1 - 2*ct2*w1w2w3 + st*theta2 - st*w1w1 + st*w1w2w3);
dRdom(6, 2) = inv_theta2*(stw3*w2w2 - stw1*w2 + ct2*theta2*w3 + cs*w1*w2 - 2*ct2*w3*w2w2);
dRdom(6, 3) = inv_theta2*(stw2*w3w3 - stw1*w3 + ct2*theta2*w2 + cs*w1*w3 - 2*ct2*w2*w3w3);
dRdom(7, 1) = inv_theta2*(stw3*w1w1 - stw1*w2 + ct2*theta2*w3 + cs*w1*w2 - 2*ct2*w3*w1w1);
dRdom(7, 2) = inv_theta2*(cs*w2w2 - 2*ct2*w1w2w3 + st*theta2 - st*w2w2 + st*w1w2w3);
dRdom(7, 3) = inv_theta2*(stw1*w3w3 - stw2*w3 + ct2*theta2*w1 + cs*w2*w3 - 2*ct2*w1*w3w3);
dRdom(8, 1) = -inv_theta2*(cs*w1w1 + 2*ct2*w1w2w3 + st*theta2 - st*w1w1 - st*w1w2w3);
dRdom(8, 2) = inv_theta2*(stw1*w2 + stw3*w2w2 + ct2*theta2*w3 - cs*w1*w2 - 2*ct2*w3*w2w2);
dRdom(8, 3) = inv_theta2*(stw1*w3 + stw2*w3w3 + ct2*theta2*w2 - cs*w1*w3 - 2*ct2*w2*w3w3);
dRdom(9, 1) = 2*ct2*inv_theta2*w1*(w1w1 + w2w2) - inv_theta2*stw1*(w1w1 + w2w2) - 2*ct2*w1;
dRdom(9, 2) = 2*ct2*inv_theta2*w2*(w1w1 + w2w2) - inv_theta2*stw2*(w1w1 + w2w2) - 2*ct2*w2;
dRdom(9, 3) = inv_theta2*w3*(w1w1 + w2w2)*(2*ct2 - st);


norm(dR1-dRdom)

%% domdR

omega = randn(3,1);
R = rodrigues(omega);

% R = rodrigues(omega/norm(omega)*pi);
% R = eye(3);
% R = R+randn(3,3)/1000;
% R = rodrigues(rodrigues(R));

R1 = R+randn(3,3)/10000;
R1 = rodrigues(rodrigues(R1));
dR = R1 - R;   

[omc,dom] = rodrigues(R) 
[om2] = rodrigues(R+dR);
om_app = omc + dom*dR(:);  
gain1 = norm(om2 - omc)/norm(om2 - om_app) 


r11=R(1,1);
r12=R(1,2);
r13=R(1,3);
r21=R(2,1);
r22=R(2,2);
r23=R(2,3);
r31=R(3,1);
r32=R(3,2);
r33=R(3,3);
theta = atan2(((r12/2 - r21/2)^2 + (r13/2 - r31/2)^2 + (r23/2 - r32/2)^2)^(1/2), r11/2 + r22/2 + r33/2 - 1/2);
xr = (r11 + r22 + r33 - 1)/2;
yrv=[r32-r23; r13-r31; r21-r12]/2;
yr = sqrt(yrv(1)^2+yrv(2)^2+yrv(3)^2);
tr2 = (r11 + r22 + r33 - 1)^2;
ant = (r12 - r21)^2 + (r13 - r31)^2 + (r23 - r32)^2;


if yr<1e-6
    if xr>0
        om=[0;0;0];
        domdR=[0 0 0 0 0 1/2 0 -1/2 0;
               0 0 -1/2 0 0 0 1/2 0 0;
               0 1/2 0 -1/2 0 0 0 0 0];
    else
        omn2=(R+eye(3))/2;
        [U,S,V]=svd(omn2);
        vec=U(:,1)/norm(U(:,1));
        theta = atan2(yr,xr);
        om = theta*vec;
        fprintf(1,'WARNING!!! Jacobian domdR (theta=pi) undefined!!!\n');
        domdR = NaN*ones(3,9);
    end
else
    vec=yrv/yr;  
    theta = atan2(yr,xr);
    om = theta*vec;
    
    domdR(1, 1) = (r23 - r32)/(ant + tr2);
    domdR(1, 2) = -((r12 - r21)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(1, 3) = -((r13 - r31)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(1, 4) = ((r12 - r21)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(1, 5) = (r23 - r32)/(ant + tr2);
    domdR(1, 6) = theta/(2*yr) - (theta*(r23 - r32)^2)/(8*yr^3) + (xr*(r23 - r32)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(1, 7) = ((r13 - r31)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(1, 8) = (theta*(r23 - r32)^2)/(8*yr^3) - theta/(2*yr) - (xr*(r23 - r32)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(1, 9) = (r23 - r32)/(ant + tr2);
    domdR(2, 1) = -(r13 - r31)/(ant + tr2);
    domdR(2, 2) = ((r12 - r21)*(r13 - r31)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(2, 3) = (theta*(r13 - r31)^2)/(8*yr^3) - theta/(2*yr) - (xr*(r13 - r31)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(2, 4) = -((r12 - r21)*(r13 - r31)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(2, 5) = -(r13 - r31)/(ant + tr2);
    domdR(2, 6) = ((r13 - r31)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(2, 7) = theta/(2*yr) - (theta*(r13 - r31)^2)/(8*yr^3) + (xr*(r13 - r31)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(2, 8) = -((r13 - r31)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(2, 9) = -(r13 - r31)/(ant + tr2);
    domdR(3, 1) = (r12 - r21)/(ant + tr2);
    domdR(3, 2) = theta/(2*yr) - (theta*(r12 - r21)^2)/(8*yr^3) + (xr*(r12 - r21)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(3, 3) = -((r12 - r21)*(r13 - r31)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(3, 4) = (theta*(r12 - r21)^2)/(8*yr^3) - theta/(2*yr) - (xr*(r12 - r21)^2)/(2*yr^2*(ant + 4*xr^2));
    domdR(3, 5) = (r12 - r21)/(ant + tr2);
    domdR(3, 6) = -((r12 - r21)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(3, 7) = ((r12 - r21)*(r13 - r31)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(3, 8) = ((r12 - r21)*(r23 - r32)*(ant*theta - 4*xr*yr + 4*theta*xr^2))/(8*yr^3*(ant + 4*xr^2));
    domdR(3, 9) = (r12 - r21)/(ant + tr2);
    
end

om
domdR
% norm(dom-domdR)
om_app = omc + domdR*dR(:);  
gain2 = norm(om2 - omc)/norm(om2 - om_app)

