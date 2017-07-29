function yy = slerp90(x,y,xx)
% SLERP90 - Spherical Linear interpolation within spans less than pi/2.
%
%   YY = SLERP90(X,Y,XX) provide NDIM(>=2) dimensional spherical
%   interpolation YY at the values of XX. X are monotonic parameters like
%   time -- N stops. Y(1:NDIM, 1:N), which will be transformed to ensure
%   every span less than pi/2, are assumed to be N unit vectors on the
%   NDIM-D sphere. XX are the positions where to interpolate. When NDIM==4,
%   this routine can be used to interpolate keys in quaternion sequence.
%
%   Algorithm: Advanced Animation And Rendering Techniques - Theory And
%   Practice, P364.
%   See also slerp, spline_interp1, rodrigues, trans_quat_axis, trans_quat_mat.

% by zpf, form BIT, 2015-7-12

bigeps = 1e-12;   % sqrt(bigeps) = 1e-6

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
h = diff(x);
assert(all(h>0) || all(h<0),'The first argument must be monotonic!');

[ndim,n]=size(y);
assert( length(x)==n && n>=2,'Number of breaks in the 1st two arguments should match and greater than 2!');
assert(ndim>=2, 'The 2nd input argument must have at least 2 rows to make interesting circle -- sphere!');
t = sqrt(sum(y.*y));
assert(all(t>0), 'Sphere vectors must not be zero vectors!');
nd1 = ones(ndim,1);
y = y./t(nd1,:);   % ensure unit vector
if h(1)<0,
    n1 = n:-1:1;
    x = x(n1);
    y = y(:,n1);
    h = -h(n1(2:end));
end;

% transform y to ensure cosom>=0
cosom = zeros(1,n-1);
p = y(:,1);
for i=2:n,
    q = y(:, i);
    idx = sum(p.*q);
    cosom(i-1) = abs(idx);
    if idx<0,
        p = -q;
        y(:, i) = p;
    else
        p = q;
    end;
end;

p = y(:,1:n-1);
q = y(:,2:n);
omega = acos(cosom);     % 0<=omega<=pi/2
sinom = sin(omega);
idx = 1-cosom>bigeps;       % 2*sin(omega/2)^2>eps, eps1<omega<=pi/2

assert(isvector(xx),'The third argument must be a vector!');
xx = xx(:)';
Np = length(xx);
[~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).

% now go to local coordinates
t = (xx-x(index))./h(index);
% local values
p = p(:, index);
q = q(:, index);
omega = omega(index);
sinom = sinom(index);
yy = zeros(ndim, Np);

id = idx(index);       % eps1<omega<=pi/2
if any(id),
    sclp = sin((1-t(id)).*omega(id))./sinom(id);
    sclq = sin(t(id).*omega(id))./sinom(id);
    yy(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

id = ~id;                % 0<=omega<=eps1
if any(id),
    sclp = 1-t(id);
    sclq = t(id);
    yy(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

return;



%% test
ndim = 3;
n = 10;
delta = 1/20;
x1= 1:n;
xx1=1:.1:n;
y1=randn(ndim,n);
ny1 = sqrt(sum(y1.*y1));
y1 = y1./ny1(ones(ndim,1),:);
p1 = y1(:,1);
for i=2:n,
    q1 = y1(:, i);
    if sum(p1.*q1)<0,
        p1= -q1;
        y1(:, i) = p1;
    else
        p1 = q1;
    end;
end;
yy1 = slerp90(x1,y1,xx1);
set(gcf,'color','w');
h = plot3(y1(1,:), y1(2,:), y1(3,:), 'ro', yy1(1,:), yy1(2,:), yy1(3,:), 'b-');
axis equal off;
axis([-1  1 -1 1 -1 1]);

% rotate(h,zdir,90)     % 旋转句柄对象

for i = 1:n,
    text(y1(1,i)+delta, y1(2,i)+delta, y1(3,i)+delta, num2str(i), 'color', 'r');
end;

filename = 'slerp90.gif';
view(3);
% [az,el]=view;
az0 = 0;
el = 20;
dt = 0.1;  % gif 的帧时间
% 沿az方向旋转一周
daz = 5;  % 每一帧方位角（azimuth）的增量
n = 360/daz;   % 旋转一周步数

for i = 1 : n-1
    view(az0 + daz*(i-1),el)
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if  i == 1,
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end;
end;
