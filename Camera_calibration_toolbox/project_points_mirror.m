function [xp,dxpdQ,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha,dxpdX] = project_points_mirror(X,Qc,Tc,hand,f,c,k,alpha)

% project_points_mirror
%
% Projects a 3D structure onto the image plane.
%
% INPUT: X: 3D structure in the world coordinate frame (3xN or 2xN matrix for N points)
%       (Qc,Tc): Rigid motion parameters between world coordinate frame and camera reference frame
%       Qc: quaternion rotation vector (4x1 vector);
%       Tc: translation vector (3x1 vector);
%       hand: handness of camera frame wrt world frame (1 or -1);
%       f: camera focal length in units of horizontal and vertical pixel units (2x1 vector)
%       c: principal point location in pixel units (2x1 vector)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%       alpha: Skew coefficient between x and y pixel (alpha = 0 <=> square pixels)
%
% OUTPUT: xp: Projected pixel coordinates (2xN matrix for N points)
%        dxpdX: Derivative of xp with respect to X(2Nx3N matrix)
%        dxpdQ: Derivative of xp with respect to Qc ((2N)x4 matrix)
%        dxpdT: Derivative of xp with respect to T ((2N)x3 matrix)
%        dxpdf: Derivative of xp with respect to f ((2N)x2 matrix if f is 2x1, or (2N)x1 matrix is f is a scalar)
%        dxpdc: Derivative of xp with respect to c ((2N)x2 matrix)
%        dxpdk: Derivative of xp with respect to k ((2N)x4 matrix)
%
% Definitions:
% Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
% The coordinate vector of P in the camera reference frame is: Xc = R*X + T
% where R is the rotation matrix corresponding to the quaternion vector Qc: R = trans_quat_mat(Qc);
% call x, y and z the 3 coordinates of Xc: x = Xc(1); y = Xc(2); z = Xc(3);
% The pinehole projection coordinates of P is [a;b] where a=x/z and b=y/z.
% call r^2 = a^2 + b^2.
% The distorted point coordinates are: xd = [xx;yy] where:
%
% xx = a * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      2*kc(3)*a*b + kc(4)*(r^2 + 2*a^2);
% yy = b * (1 + kc(1)*r^2 + kc(2)*r^4 + kc(5)*r^6)      +      kc(3)*(r^2 + 2*b^2) + 2*kc(4)*a*b;
%
% The left terms correspond to radial distortion (6th degree), the right terms correspond to tangential distortion
%
% Finally, convertion into pixel coordinates: The final pixel coordinates vector xp=[xxp;yyp] where:
%
% xxp = f(1)*(xx + alpha*yy) + c(1)
% yyp = f(2)*yy + c(2)
%
%
% Important function related with the program:
% trans_quat_mat: transform rotation matrix to quternion vector, or vice versa.
% dAB: derivatives of AB wrt A, and AB wrt B. 
% rigid_trans: Computes the rigid transformation of a given structure
%  See also project_points2.


if nargin <8,
    alpha = 0;
    if nargin < 7,
        k = zeros(5,1);
        if nargin < 6,
            c = zeros(2,1);
            if nargin < 5,
                f = ones(2,1);
                if nargin < 4,
                    hand = 1;
                    if nargin < 3,
                        Tc = zeros(3,1);
                        if nargin < 2,
                            Qc = [zeros(3,1); 1];
                            if nargin < 1,
                                error('Need at least a 3D structure to project!');
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

[m,n] = size(X);

if nargout > 1,
    [Y,dYdQ,dYdT,dYdX] = rigid_trans(X,Qc,Tc,hand);
else
    Y = rigid_trans(X,Qc,Tc,hand);
end;

inv_Z = 1./Y(3,:);
x = Y(1:2,:) .*inv_Z (ones(2,1), :) ;
% Add distortion:
r2 = x(1,:).^2 + x(2,:).^2;
r4 = r2.^2;
r6 = r2.^3;
% Radial distortion:
cdist = 1 + k(1) * r2 + k(2) * r4 + k(5) * r6;
xd1 = x .* (ones(2,1)*cdist);
% tangential distortion:
a1 = 2.*x(1,:).*x(2,:);
a2 = r2 + 2*x(1,:).^2;
a3 = r2 + 2*x(2,:).^2;

delta = [k(3)*a1 + k(4)*a2 ;
    k(3) * a3 + k(4)*a1];

xd2 = xd1 + delta;
% Add Skew:
xd3 = [xd2(1,:)+alpha*xd2(2,:); xd2(2,:)];
% Pixel coordinates:
if length(f)>1,
    xp = xd3 .* (f(:) * ones(1,n))  +  c(:)*ones(1,n);
else
    xp = f * xd3 + c*ones(1,n);
end;

if nargout<2,
    return;
end;

n2 = n*2;
% dxdX, dxdQ, dxdT
xn = (-x(1,:) .* inv_Z)';
yn = (-x(2,:) .* inv_Z)';
dxdQ = zeros(n2,4);
dxdT = zeros(n2,3);
dxdQ(1:2:end,:) = (inv_Z'*ones(1,4)) .* dYdQ(1:3:end,:) + xn(:, ones(1,4)) .* dYdQ(3:3:end,:);
dxdQ(2:2:end,:) = (inv_Z'*ones(1,4)) .* dYdQ(2:3:end,:) + yn(:, ones(1,4)) .* dYdQ(3:3:end,:);

dxdT(1:2:end,:) = (inv_Z'*ones(1,3)) .* dYdT(1:3:end,:) + xn(:, ones(1,3)) .* dYdT(3:3:end,:);
dxdT(2:2:end,:) = (inv_Z'*ones(1,3)) .* dYdT(2:3:end,:) + yn(:, ones(1,3)) .* dYdT(3:3:end,:);


% dr2dQ, dr2dT
dr2dQ = 2*(x(1,:)'*ones(1,4)).* dxdQ(1:2:end,:) + 2*(x(2,:)'*ones(1,4)) .* dxdQ(2:2:end,:);
dr2dT = 2*(x(1,:)'*ones(1,3)).* dxdT(1:2:end,:) + 2*(x(2,:)'*ones(1,3)) .* dxdT(2:2:end,:);

% dcdistdQ, dcdistdT, dcdistdk
dcdistdQ = k(1) * dr2dQ + ((k(2)*2*r2' + k(5)*3*r4')*ones(1,4)).* dr2dQ;
dcdistdT = k(1) * dr2dT + ((k(2)*2*r2' + k(5)*3*r4')*ones(1,3)).*dr2dT;
dcdistdk = [ r2' r4' zeros(n,2) r6'];

% dxd1dQ, dxd1dT, dxd1dk
dxd1dQ = zeros(n2,4);
dxd1dT = zeros(n2,3);
dxd1dk = zeros(n2,5);
dxd1dQ(1:2:end,:) = (x(1,:)'*ones(1,4)) .* dcdistdQ;
dxd1dQ(2:2:end,:) = (x(2,:)'*ones(1,4)) .* dcdistdQ;
coeff = reshape([cdist;cdist],n2,1);
dxd1dQ = dxd1dQ + coeff(:,ones(1,4)).* dxdQ;

dxd1dT(1:2:end,:) = (x(1,:)'*ones(1,3)) .* dcdistdT;
dxd1dT(2:2:end,:) = (x(2,:)'*ones(1,3)) .* dcdistdT;
dxd1dT = dxd1dT + coeff(:,ones(1,3)).* dxdT;

dxd1dk(1:2:end,:) = (x(1,:)'*ones(1,5)) .* dcdistdk;
dxd1dk(2:2:end,:) = (x(2,:)'*ones(1,5)) .* dcdistdk;

% ddeltadQ, ddeltadT, ddeltadk
%ddeltadx = zeros(2*n,2*n);
aa = (2*k(3)*x(2,:)+6*k(4)*x(1,:))';
bb = (2*k(3)*x(1,:)+2*k(4)*x(2,:))';
cc = (6*k(3)*x(2,:)+2*k(4)*x(1,:))';
ddeltadQ = zeros(n2,4);
ddeltadQ(1:2:end,:) = aa(:,ones(1,4)) .* dxdQ(1:2:end,:) + bb(:,ones(1,4)) .* dxdQ(2:2:end,:);
ddeltadQ(2:2:end,:) = bb(:,ones(1,4)) .* dxdQ(1:2:end,:) + cc(:,ones(1,4)) .* dxdQ(2:2:end,:);

ddeltadT = zeros(n2,3);
ddeltadT(1:2:end,:) = aa(:,ones(1,3)) .* dxdT(1:2:end,:) + bb(:,ones(1,3)) .* dxdT(2:2:end,:);
ddeltadT(2:2:end,:) = bb(:,ones(1,3)) .* dxdT(1:2:end,:) + cc(:,ones(1,3)) .* dxdT(2:2:end,:);

ddeltadk = zeros(n2,5);
ddeltadk(1:2:end,3) = a1';
ddeltadk(1:2:end,4) = a2';
ddeltadk(2:2:end,3) = a3';
ddeltadk(2:2:end,4) = a1';

% dxd2dQ, dxd2dT, dxd2dk
dxd3dQ = dxd1dQ + ddeltadQ ;
dxd3dT = dxd1dT + ddeltadT;
dxd3dk = dxd1dk + ddeltadk ;

% Compute: dxd3dQ, dxd3dT, dxd3dk, dxd3dalpha
dxd3dQ(1:2:n2,:) = dxd3dQ(1:2:n2,:)  + alpha*dxd3dQ(2:2:n2,:);
dxd3dT(1:2:n2,:) = dxd3dT(1:2:n2,:) + alpha*dxd3dT(2:2:n2,:);
dxd3dk(1:2:n2,:) = dxd3dk(1:2:n2,:) + alpha*dxd3dk(2:2:n2,:);
dxd3dalpha = zeros(n2,1);
dxd3dalpha(1:2:n2,:) = xd2(2,:)';

% dxpdQ, dxpdT, dxpdk, dxpdalpha, dxpdf, dxpdc
if length(f)>1,
    ff = reshape(f(:)*ones(1,n),n2,1);
    dxpdQ = ff(:,ones(1,4)) .* dxd3dQ;
    dxpdT = ff(:,ones(1,3)) .* dxd3dT;
    dxpdk = ff(:,ones(1,5)) .* dxd3dk;
    dxpdalpha = ff .* dxd3dalpha;
    dxpdf = zeros(n2,2);
    dxpdf(1:2:end,1) = xd3(1,:)';
    dxpdf(2:2:end,2) = xd3(2,:)';
else
    dxpdQ = f  * dxd3dQ;
    dxpdT = f * dxd3dT;
    dxpdk = f  * dxd3dk;
    dxpdalpha = f .* dxd3dalpha;
    dxpdf = xd3(:);
end;

dxpdc = zeros(n2,2);
dxpdc(1:2:end,1) = ones(n,1);
dxpdc(2:2:end,2) = ones(n,1);

if nargout<8,
    return;
end;

nm = n*m;
if n>5,
    dxdX = sparse([],[],[],n2,nm,n*4);
else
    dxdX = zeros(n2,nm); 
end;
dxd1dX = dxdX;
ddelta_xdX =dxdX;
dxdX(1:2:end,:) = (inv_Z'*ones(1,nm)) .* dYdX(1:3:end,:) + xn(:, ones(1,nm)) .* dYdX(3:3:end,:);
dxdX(2:2:end,:) = (inv_Z'*ones(1,nm)) .* dYdX(2:3:end,:) + yn(:, ones(1,nm)) .* dYdX(3:3:end,:);
dr2dX = 2*(x(1,:)'*ones(1,nm)) .* dxdX(1:2:end,:) + 2*(x(2,:)'*ones(1,nm)) .* dxdX(2:2:end,:);
dcdistdX = k(1) * dr2dX + ((k(2)*2*r2' + k(5)*3*r4')*ones(1,nm)).* dr2dX;
dxd1dX(1:2:end,:) = (x(1,:)'*ones(1,nm)) .* dcdistdX;
dxd1dX(2:2:end,:) = (x(2,:)'*ones(1,nm)) .* dcdistdX;
dxd1dX = dxd1dX + coeff(:,ones(1,nm)).* dxdX;
ddelta_xdX(1:2:end,:) = aa(:,ones(1,nm)) .* dxdX(1:2:end,:) + bb(:,ones(1,nm)) .* dxdX(2:2:end,:);
ddelta_xdX(2:2:end,:) = bb(:,ones(1,nm)) .* dxdX(1:2:end,:) + cc(:,ones(1,nm)) .* dxdX(2:2:end,:);
dxd3dX = dxd1dX + ddelta_xdX ;
dxd3dX(1:2:n2,:) = dxd3dX(1:2:n2,:)  + alpha*dxd3dX(2:2:n2,:);
if length(f)>1,
    ff = reshape(f(:)*ones(1,n),n2,1);
    dxpdX = ff(:,ones(1, nm)) .* dxd3dX;
else
    dxpdX = f  * dxd3dX;
end;

return;



%% Test of the Jacobians:
m = 2;
n = 1000;
X = 10*randn(m,n);
om = randn(3,1);
q = trans_quat_axis(om);
T = [10*randn(2,1);40];
hand = sign(randn(1));
f = 1000*rand(2,1);
c = 1000*randn(2,1);
k = 0.5*randn(5,1);
alpha = 0.01*randn(1,1);
[x,dxdq,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdX] = project_points_mirror(X,q,T,hand,f,c,k,alpha);
dX = randn(m,n)/1000;
dom = randn(3,1)/1000;
dq = trans_quat_axis(om+dom)-q;
dT = norm(T)*randn(3,1)/1000;
df = norm(f)*randn(2,1)/1000;
dc = norm(c)*randn(2,1)/100;
dk = norm(k)*randn(5,1)/1000;
dalpha = norm(k)*randn(1,1)/1000;
x2 = project_points_mirror(X+dX,q+dq,T+dT,hand,f+df,c+dc,k+dk,alpha+dalpha);
x_pred = x +reshape(dxdX*dX(:),2,n)+ reshape(dxdq*dq,2,n)+reshape(dxdT*dT,2,n)+reshape(dxdf*df,2,n)+reshape(dxdc*dc,2,n)+reshape(dxdk*dk,2,n)+reshape(dxdalpha*dalpha,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on X
x2 = project_points_mirror(X+dX,q,T,hand,f,c,k,alpha);
x_pred = x + reshape(dxdX * dX(:),2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on q
x2 = project_points_mirror(X,q+dq,T,hand,f,c,k,alpha);
x_pred = x + reshape(dxdq * dq,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on T
x2 = project_points_mirror(X,q,T+dT,hand,f,c,k,alpha);
x_pred = x + reshape(dxdT * dT,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on f
[x2] = project_points_mirror(X,q,T,hand,f+df,c,k,alpha);
x_pred = x + reshape(dxdf * df,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on c
[x2] = project_points_mirror(X,q,T,hand,f,c+dc,k,alpha);
x_pred = x + reshape(dxdc * dc,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)


%% Test on k
[x2] = project_points_mirror(X,q,T,hand,f,c,k+dk,alpha);
x_pred = x + reshape(dxdk * dk,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)

%% Test on alpha
[x2] = project_points_mirror(X,q,T,hand,f,c,k,alpha+dalpha);
x_pred = x + reshape(dxdalpha * dalpha,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)
