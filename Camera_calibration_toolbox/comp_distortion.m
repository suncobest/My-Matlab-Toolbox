function x_comp  = comp_distortion(x_dist,k)

%       [x_comp] = comp_distortion(x_dist,k);
%
%       compensates the radial distortion of the camera
%       on the image plane. Model From Oulu university.
%       For more information about the distortion model,
%       check the forward projection mapping function: project_points.m
%
%       x_dist : Distorted point coordinates in the image plane (2xN matrix)
%       k: Distortion coefficients (radial and tangential
%       x_comp : The image plane points after correction for the distortion (2xN matrix)
%
%       NOTE : This compensation has to be done after the substraction
%              of the center of projection, and division by the focal length.
%
%       If length(k) == 1:
%       Solve for cubic roots using method from Numerical Recipes in C++ 3rd Ed.
%       pages 228-229.
%
%       rd=ru*(1+k*ru^2) --- x^3+(1/k)*x-rd/k=0:  x=ru, a=0,b=1/k,c=-rd/k

assert(ismatrix(x_dist) && size(x_dist,1)==2,'The dimension of the points should be 2xN!');
m = length(k);
if m < 5,
    k = [k(:); zeros(5-m,1)];
end;

if norm(k)==0,
    x_comp = x_dist;
    return;
elseif norm(k(2:end))==0,
    k1 = k(1);
    ph = atan2(x_dist(2,:),x_dist(1,:));
    Q = -1/(3*k1);
    R = -x_dist(1,:)./(2*k1*cos(ph));          % R=-rd/(2*k1), rd=x/cos(ph)
    R2 = R.^2;
    Q3 = Q^3;

    if k1 < 0,   % R2<Q3
        % this works in all practical situations (it starts failing for very large
        % values of k1)
        th = acos(R./sqrt(Q3));                    % in test: 0<th<pi/2
        r = -2*sqrt(Q)*cos((th-2*pi)/3);           % -2*pi/3<(th-2*pi)/3<-pi/2, r>0
    else
        % note: this always works, even for ridiculous values of k1
        A = (sqrt(R2-Q3)-R).^(1/3);
        B = Q./A;
        r = (A+B);
    end;
    x_comp = real([r.*cos(ph); r.*sin(ph)]);
    return;
end;

k1 = k(1);
k2 = k(2);
k3 = k(5);
p1 = k(3);
p2 = k(4);

x_comp = x_dist;   % initial guess
for kk=1:10,
    r_2 = sum(x_comp.^2,1);
    k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
    delta_x = [2*p1*x_comp(1,:).*x_comp(2,:) + p2*(r_2 + 2*x_comp(1,:).^2);
    p1 * (r_2 + 2*x_comp(2,:).^2)+2*p2*x_comp(1,:).*x_comp(2,:)];
    x_comp = (x_dist - delta_x)./(ones(2,1)*k_radial);
end;

return;


%% test
n = 10;
lk = randi(5);
kc = randn(lk,1)/10;
xc = randn(2,n);
xd = apply_distortion(xc,kc);
xc1 = comp_distortion(xd,kc);
norm(xc1-xc)
