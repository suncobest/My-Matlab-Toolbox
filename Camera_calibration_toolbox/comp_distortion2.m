function [x_comp]  = comp_distortion2(x_dist,k)

%       [x_comp] = comp_distortion(x_dist,k);
%
%       compensates the radial distortion of the camera
%       on the image plane.
%
%       x_dist : the image points got without considering the
%                radial distortion.
%       k: Radial distortion factor
%
%       x : The image plane points after correction for the distortion
%
%       x and x_dist are 2xN arrays
%
%       NOTE : This compensation has to be done after the substraction
%              of the center of projection, and division by the focal
%              length.
%
%  Solve for cubic roots using method from Numerical Recipes in C 2nd Ed.
%  pages 184-185.
%
% rd=ru*(1+k*ru^2) --- x^3+(1/k)*x-rd/k=0:  x=ru, a=0,b=1/k,c=-rd/k

% California Institute of Technology
% (c) Jean-Yves Bouguet - April 27th, 1998
% fully checked! JYB, april 27th, 1998 - 2am

if k == 0,
    x_comp = real(x_dist);
    return;
end;

[two,N] = size(x_dist);
if (two ~= 2 ),
    error('The dimension of the points should be 2xN!');
end;
ph = atan2(x_dist(2,:),x_dist(1,:));
Q = -1/(3*k);
R = -x_dist(1,:)./(2*k*cos(ph));          % R=-rd/(2*k), rd=x/cos(ph)
R2 = R.^2;
Q3 = Q^3;

if k < 0,
   % this works in all practical situations (it starts failing for very large
   % values of k)
   th = acos(R./sqrt(Q3));                    % in test: 0<th<pi/2
   r = -2*sqrt(Q)*cos((th-2*pi)/3);           % -2*pi/3<(th-2*pi)/3<-pi/2, r>0
else
   % note: this always works, even for ridiculous values of k
   A = (sqrt(R2-Q3)-R).^(1/3);
   B = Q*(1./A);
   r = (A+B);
end;
x_comp = real([r.*cos(ph); r.*sin(ph)]);
