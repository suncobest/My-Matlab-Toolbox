function [Y,dYdX] = spherical2cartesian(X)

% SPHERICAL2CARTESIAN compute the corresponding Cartesian coordinates from the
% original spherical coordinates, where X=[r;theta;phi], Y=[x;y;z]. Theta is the
% angle from z axis to the vector, and phi is the angle from x axis to the vector.
%
% Given X=[r;theta;phi], and Y=[x;y;z], where
%
%  x = r*sin(theta)*cos(phi);   y = r*sin(theta)*sin(phi);   z = r*cos(theta);
%
%  The Jacobian matrix of dYdX is:
%
%         --                                                                      --
%         | sin(theta)*cos(phi)    r*cos(theta)*cos(phi)    -r*sin(theta)*sin(phi) |
%         |                                                                        |
%  dYdX = | sin(theta)*sin(phi)    r*cos(theta)*sin(phi)     r*sin(theta)*cos(phi) |
%         |                                                                        |
%         |    cos(theta)              -r*sin(theta)                     0         |
%         --                                                                      --
%
%  See also cartesian2spherical, cartesian2cylindrical, cylindrical2cartesian.

%  By ZPF @ZVR, 2017-7-20


assert(ismatrix(X) && size(X,1)==3, 'Unexpected dimension of input matrix!');

r = X(1,:);
assert(all(r>=0), 'Radius are assumed to be nonnegative!');
stheta = sin(X(2,:));
ctheta = cos(X(2,:));
sphi = sin(X(3,:));
cphi = cos(X(3,:));

sthetacphi = stheta.*cphi;
sthetasphi = stheta.*sphi;

x = r.*sthetacphi;
y = r.*sthetasphi;
z = r.*ctheta;
Y = [x; y; z];

if nargout > 1,
    zcphi = z.*cphi;
    zsphi = z.*sphi;
    rstheta = r.*stheta;

    n = size(X,2);
    n3 = n*3;
    if n>5;
        dYdX = sparse([],[],[],n3,n3,n3*3);
    else
        dYdX = zeros(n3);
    end;
    for i=1:n,
        j = (i-1)*3;
        dYdX(j+1:j+3,j+1:j+3) = [sthetacphi(i), zcphi(i), -y(i);  sthetasphi(i), zsphi(i), x(i); ctheta(i), -rstheta(i), 0];
    end;
end;

return;



%% Test of the Jacobians:
m = 1000;
phi = randn(1,m);
phi = mod(phi+pi,2*pi)-pi;
XX = [rand(1,m)*5+0.5; rand(1,m)*pi; phi];
norm(XX-cartesian2spherical(spherical2cartesian(XX)))

[YY,dYYdXX] = spherical2cartesian(XX);
dXX = randn(3,m)/100;
YY1 = spherical2cartesian(XX+dXX);
YY_pred = YY + reshape(dYYdXX*dXX(:), 3, m);
gain = norm(YY1-YY)/norm(YY1-YY_pred)
