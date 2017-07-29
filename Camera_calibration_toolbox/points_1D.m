function [X,dXdX0,dXdthetaphi,dXdx] = points_1D(X0,thetaphi,lamda)

% POINTS_1D compute one dimensional points with origin at X0, direction angle of
% thetaphi=[theta; phi], and 1D coordinates of lamda=[lamda1,lamda2,...,lamdan].
% [theta; phi] isthe same as the direction of Spherical coordinate system. Theta is the
% angle from Z axis to the vector, and phi is the angle from X axis to the vector;
%
%  X = X0+r*lamda;
%
%  See also points_1D2, cartesian2spherical, project_points_mirror2, rigid_refmotion.

% By ZPF @ZVR, 2017-7-19


lamda = lamda(:)';
n = length(lamda);

theta = thetaphi(1);
phi = thetaphi(2);

ctheta = cos(theta);
stheta = sin(theta);
cphi = cos(phi);
sphi = sin(phi);

% direction of the line.
r = [stheta*cphi; stheta*sphi; ctheta];

X = X0(:,ones(1,n))+r*lamda;

if nargout > 1,
    dXdX0 = repmat(eye(3),n,1);
    if nargout > 2,
        n3 = 3*n;
        dXdr = zeros(n3,3);
        for i=1:n,
            j = (i-1)*3;
            dXdr(j+1:j+3,:) = lamda(i)*eye(3);
        end;
        drdthetaphi = [ctheta*cphi,  -stheta*sphi;
                       ctheta*sphi,  stheta*cphi;
                         -stheta,         0];
        dXdthetaphi = dXdr*drdthetaphi;
        if nargout > 3,
            dXdx = repmat([r; zeros(n3,1)],1,n);
            dXdx = reshape(dXdx(1:n3*n),n3,n);
        end;
    end;
end;

return;



%% Test of the Jacobians:
m = 10;
xx = randn(m,1);
tp = randn(2,1);
tp(1) = 0;
XX0 = 10*randn(3,1);
[XX,dXXdXX0,dXXdtp,dXXdxx] = points_1D(XX0,tp,xx);
dXX0 = randn(3,1)/100;
dtp = randn(2,1)/1000;
dxx = randn(m,1)/100;
XX1 = points_1D(XX0+dXX0,tp+dtp,xx+dxx);
XX_pred = XX + reshape(dXXdXX0*dXX0 + dXXdtp*dtp + dXXdxx*dxx, 3, m);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on XX0
XX1 = points_1D(XX0+dXX0,tp,xx);
XX_pred = XX + reshape(dXXdXX0*dXX0, 3, m);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on tp
XX1 = points_1D(XX0,tp+dtp,xx);
XX_pred = XX + reshape(dXXdtp*dtp, 3, m);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on xx
XX1 = points_1D(XX0,tp,xx+dxx);
XX_pred = XX + reshape(dXXdxx*dxx, 3, m);
gain = norm(XX1-XX)/norm(XX1-XX_pred)
