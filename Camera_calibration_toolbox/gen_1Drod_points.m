function [X,dXdXori,dXdthphi,dXdlen] = gen_1Drod_points(Xori,thetaphi,rodlen)

% GEN_1DROD_POINTS compute one dimensional points with origin at Xori, direction angle
% of thetaphi=[theta; phi], and 1D coordinates of rodlen=[x1,x2,...,xn].
% [theta; phi] isthe same as the direction of Spherical coordinate system. Theta is the
% angle from Z axis to the vector, and phi is the angle from X axis to the vector;
%
%  X = Xori+Xn*rodlen;
%
%  See also points_1D2, cartesian2spherical, project_points_mirror2, rigid_refmotion.

% By ZPF @ZVR, 2017-7-19


[m,N] = size(Xori);
assert(ismatrix(Xori) && m==3,'Unexpected dimension of the 1st input!');
[m,n] = size(thetaphi);
assert(ismatrix(thetaphi) && m==2,'Unexpected dimension of the 2nd input!');
if N ~= n,
    if N==1,
        Xori = Xori(:,ones(1,n));
        N = n;
    elseif n==1,
        thetaphi = thetaphi(:,ones(1,N));
    else
        error('The columns of the 1st two variables do not match!');
    end;
end;

rodlen = rodlen(:)';
np1D = length(rodlen);

ctheta = cos(thetaphi(1,:));
stheta = sin(thetaphi(1,:));
cphi = cos(thetaphi(2,:));
sphi = sin(thetaphi(2,:));

% direction of the line.
Xn = [stheta.*cphi; stheta.*sphi; ctheta];

X = reshape(permute(Xori+Xn.*reshape(rodlen,[1,1,np1D]), [1,3,2]), 3,[]);




if nargout > 1,
    dXdXori = repmat(eye(3),np1D,1);
    if nargout > 2,
        n3 = 3*np1D;
        dXdr = zeros(n3,3);
        for i=1:np1D,
            j = (i-1)*3;
            dXdr(j+1:j+3,:) = rodlen(i)*eye(3);
        end;
        dXndthphi = [ctheta*cphi,  -stheta*sphi;
                       ctheta*sphi,  stheta*cphi;
                         -stheta,         0];
        dXdthphi = dXdr*dXndthphi;
        if nargout > 3,
            dXdlen = repmat([r; zeros(n3,1)],1,np1D);
            dXdlen = reshape(dXdlen(1:n3*np1D),n3,np1D);
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
