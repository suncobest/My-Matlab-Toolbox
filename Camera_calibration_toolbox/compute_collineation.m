function [H,Hnorm,inv_Hnorm] = compute_collineation (a00, a10, a11, a01);

% new formalism using homographies

a00 = a00 / a00(3);
a10 = a10 / a10(3);
a11 = a11 / a11(3);
a01 = a01 / a01(3);


% Prenormalization of point coordinates (very important):
% (Affine normalization)

ax = [a00(1);a10(1);a11(1);a01(1)];
ay = [a00(2);a10(2);a11(2);a01(2)];

mxx = mean(ax);
myy = mean(ay);
ax = ax - mxx;
ay = ay - myy;

scxx = mean(abs(ax));
scyy = mean(abs(ay));


Hnorm = [1/scxx 0 -mxx/scxx;0 1/scyy -myy/scyy;0 0 1];
inv_Hnorm = [scxx 0 mxx ; 0 scyy myy; 0 0 1];


a00n = Hnorm*a00;
a10n = Hnorm*a10;
a11n = Hnorm*a11;
a01n = Hnorm*a01;


% Computation of the vanishing points:

V1n = cross(cross(a00n,a10n),cross(a01n,a11n));
V2n = cross(cross(a00n,a01n),cross(a10n,a11n));


% Normalizaion of the vanishing points:

V1n = V1n/norm(V1n);
V2n = V2n/norm(V2n);


% Closed-form solution of the coefficients:  a10n=a00n+alpha_x*V1n; a01n=a00n+alpha_y*V2n; 
% cross product is more accurate than dot product --- alpha_x=V1n'*(a10n-a00n), alpha_y=V2n'*(a10n-a00n); 
% three points on the same line(a00n,a10n,V1n) --- cross(a00n,a10n)=alpha_x*cross(a10n,V1n)

alpha_x = (a10n(2)*a00n(1) - a10n(1)*a00n(2))/(V1n(2)*a10n(1)-V1n(1)*a10n(2));      % cross(a00n,a10n)+alpha_x*cross(V1n,a10n)=0

alpha_y = (a01n(2)*a00n(1) - a01n(1)*a00n(2))/(V2n(2)*a01n(1)-V2n(1)*a01n(2));      % cross(a00n,a01n)+alpha_y*cross(V2n,a01n)=0


% Remaining Homography: affinity (6DOF --- can be solved by three pair of points)

Hrem = [alpha_x*V1n  alpha_y*V2n a00n];         % [a10n,a01n,a00n] = Hrem*[1 0 0; 0 1 0; 1 1 1];


% Final homography:

H = inv_Hnorm*Hrem;

