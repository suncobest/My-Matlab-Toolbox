function [Y,dYdX] = cartesian2spherical(X)

% CARTESIAN2SPHERICAL compute the corresponding spherical coordinates from the
% original Cartesian coordinates, where X=[x;y;z], Y=[r;theta;phi]. Theta is the
% angle from z axis to the vector, and phi is the angle from x axis to the vector.
%
% Given X=[x;y;z], and Y=[r;theta;phi], then
%
%  r = sqrt(x^2+y^2+z^2);   theta = arccos(z/r);    phi = arctan(y/x);
%
%  The Jacobian matrix of dYdX is:
%
%         --                                                            --
%         |       x/r                  y/r                    z/r        |
%         |                                                              |
%  dYdX = | xz/[(x^2+y^2)r^2]    yz/[(x^2+y^2)r^2]    -sqrt(x^2+y^2)/r^2 |
%         |                                                              |
%         |  -y/(x^2+y^2)           x/(x^2+y^2)                0         |
%         --                                                            --
%
%  See also spherical2cartesian, cartesian2cylindrical, cylindrical2cartesian.

%  By ZPF @ZVR, 2017-7-20


assert(ismatrix(X) && size(X,1)==3, 'Unexpected dimension of input matrix!');

X2 = X.*X;
r2 = sum(X2,1);
r = sqrt(r2);

% theta is in [0,pi]
ctheta = X(3,:)./r;
theta = acos(ctheta);

% -pi<=phi<=pi; through mod(a+pi,2*pi)-pi, a can be shift a in [-pi,pi].
phi = atan2(X(2,:), X(1,:));
Y = [r; theta; phi];

if nargout > 1,
    drdxyz = X./r(ones(3,1),:);
    x2py2 = r2-X2(3,:);
    sqx2py2 = sqrt(x2py2);
    dthetadx = X(1,:).*X(3,:)./(r2.*sqx2py2);
    dthetady = X(2,:).*X(3,:)./(r2.*sqx2py2);
    dthetadz = -sqx2py2./r2;
    dphidx = -X(2,:)./x2py2;
    dphidy = X(1,:)./x2py2;
    % dphidz = 0;

    n = size(X,2);
    n3 = n*3;
    if n>5;
        dYdX = sparse([],[],[],n3,n3,n3*3);
    else
        dYdX = zeros(n3);
    end;
    for i=1:n,
        j = (i-1)*3;
        dYdX(j+1:j+3,j+1:j+3) = [drdxyz(:,i)'; dthetadx(i), dthetady(i), dthetadz(i); dphidx(i), dphidy(i), 0];
    end;
end;

return;



%% Test of the Jacobians:
m = 1000;
XX = randn(3,m);
norm(XX-spherical2cartesian(cartesian2spherical(XX)))

[YY,dYYdXX] = cartesian2spherical(XX);
dXX = randn(3,m)/100;
YY1 = cartesian2spherical(XX+dXX);
YY_pred = YY + reshape(dYYdXX*dXX(:), 3, m);
dYY = YY1-YY;
dYYp = YY1-YY_pred;
id = dYY(3,:)>pi;
dYY(3,id) = dYY(3,id)-2*pi;
id = dYY(3,:)<-pi;
dYY(3,id) = dYY(3,id)+2*pi;

id = dYYp(3,:)>pi;
dYYp(3,id) = dYYp(3,id)-2*pi;
id = dYYp(3,:)<-pi;
dYYp(3,id) = dYYp(3,id)+2*pi;
gain = norm(dYY)/norm(dYYp)

