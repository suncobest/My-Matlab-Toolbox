function yy = squad(x,y,xx)
% SQUAD - Spherical Spline Quaternion interpolation.
%
%   YY = SQUAD(X,Y,XX) provide spherical spline quaternion interpolation YY
%   at the values of XX. X are monotonic parameters like time -- N
%   stops. Y(1:4, 1:N), which will be transformed to ensure every successive
%   span angle less than pi/2, are assumed to be N unit vectors on the 4D
%   sphere. XX are the positions where to interpolate.
%
%   Algorithm: Quaternions, Interpolation and Animation, p51.
%   See also slerp90, slerp, spline_interp1, rodrigues, trans_quat_axis, trans_quat_mat.

% by zpf, form BIT, 2015-7-22

bigeps = 1e-12;   % sqrt(bigeps) = 1e-6

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
h = diff(x);
assert(all(h>0) || all(h<0),'The first argument must be monotonic!');

[ndim,n]=size(y);
assert( length(x)==n && n>=2,'Number of breaks in the 1st two arguments should match and greater than 2!');
assert(ndim==4, 'The 2nd input argument must have 4 rows (quaternions)!');
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
omega = acos(cosom);     % 0<=omega<=pi/2
sinom = sin(omega);
idx = 1-cosom>bigeps;       % 2*sin(omega/2)^2>eps, eps1<omega<=pi/2

% caculate the auxiliary points s
p = y(:, 1:n-2);
s = y(:, 2:n-1);
q = y(:, 3:n);
invs = [-s(1:3, :); s(4, :)];
s = quatmul(s, quatexp(-(quatlog(quatmul(invs,p,1))+quatlog(quatmul(invs,q,1)))/4),1);
s = [y(:,1), s, y(:,n)];

assert(isvector(xx),'The third argument must be a vector!');
xx = xx(:)';
Np = length(xx);
[~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).
% now go to local coordinates
t = (xx-x(index))./h(index);

% % p=slerp(p1,p2,t)--p = slerp(x, y, xx);
p = y(:,1:n-1);
q = y(:,2:n);
p = p(:, index);
q = q(:, index);
omega = omega(index);
sinom = sinom(index);
y = zeros(ndim, Np);

id = idx(index);       % eps1<omega<=pi/2
if any(id),
    sclp = sin((1-t(id)).*omega(id))./sinom(id);
    sclq = sin(t(id).*omega(id))./sinom(id);
    y(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

id = ~id;                % 0<=omega<=eps1
if any(id),
    sclp = 1-t(id);
    sclq = t(id);
    y(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

% % q=slerp(s1,s2,t)--q = slerp(x, s, xx);
p = s(:,1:n-1);
q = s(:,2:n);
cosom = sum(p.*q);               % 0<=omega<=pi, 0=<omega/2<=pi/2
id1 = 1+cosom>bigeps;       % 2*cos(omega/2)^2 > eps, 0=<omega/2<pi/2-eps1
id2 = 1-cosom>bigeps;       % 2*sin(omega/2)^2 > eps, eps1<omega/2<=pi/2
omega = acos(cosom);
sinom = sin(omega);

pl = p(:, index);
ql = q(:, index);
omega = omega(index);
sinom = sinom(index);
yy = zeros(ndim, Np);

id = id1 & id2;          % 2*eps1<omega<pi-2*eps1
id = id(index);
if any(id),
    sclp = sin((1-t(id)).*omega(id))./sinom(id);
    sclq = sin(t(id).*omega(id))./sinom(id);
    yy(:,id) = sclp(nd1,:).*pl(:,id)+sclq(nd1,:).*ql(:,id);
end;

id = ~id2(index);       % 0<=omega<=2*eps1
if any(id),
    sclp = 1-t(id);
    sclq = t(id);
    yy(:,id) = sclp(nd1,:).*pl(:,id)+sclq(nd1,:).*ql(:,id);
end;

id = ~id1(index);         % pi-2*eps1<=omega<=pi, two antipole points span half the equator
if any(id),
    idx = ~id1;
    flag = 1;
    xx = 0;
    bigeps = 1e-4;
    while flag,                % make sure y is far enough away from both p and q
        while ~all(xx),      % make sure y~=0
            y = randn(ndim,sum(idx));
            xx = sqrt(sum(y.*y));
        end;
        y = y./xx(nd1,:);   % ensure unit vector
        xx = 0;
        om1 = acos(sum(p(:, idx).*y));
        om2 = acos(sum(q(:, idx).*y));
        om = om1+om2;
        assert(all(abs(om-pi)<bigeps), 'Antipoles should span pi!');
        flag = ~all(om1>bigeps & om1<pi-bigeps);
    end;
    yyy = zeros(ndim,n-1);
    yyy(:, idx) = y;
    yyy = yyy(:, index);
    y = yyy(:,id);
    
    omega = zeros(1,n-1); 
    omega(idx) = om1;
    omega = omega(index);
    om1 = omega(id);
    
    omega = zeros(1,n-1);
    omega(idx) = om2;
    omega = omega(index);
    om2 = omega(id);
    om = om1+om2;
    
    pl = pl(:,id);
    ql = ql(:,id);
    tom1 = om1./om; 
    tom2 = 1-tom1;
    ti = t(id);
    yyy = zeros(ndim, sum(id));
    
    idx = ti<tom1;
    ty = ti(idx)./tom1(idx);
    sinom = sin(om1);
    sclp = sin((1-ty).*om1(idx))./sinom(idx);
    sclq = sin(ty.*om1(idx))./sinom(idx);
    yyy(:,idx) = sclp(nd1,:).*pl(:,idx) + sclq(nd1,:).*y(:,idx);
    
    idx = ~idx;
    ty = (ti(idx)-tom1(idx))./tom2(idx);
    sinom = sin(om2);
    sclp = sin((1-ty).*om2(idx))./sinom(idx);
    sclq = sin(ty.*om2(idx))./sinom(idx);
    yyy(:,idx) = sclp(nd1,:).*y(:,idx) + sclq(nd1,:).*ql(:,idx);
    
    yy(:,id) = yyy;
end;

% % yy = slerp(p, q, 2t(1-t));
p = y;
q = yy;
cosom = sum(p.*q);
omega = acos(cosom);     % 0<=omega<=pi/2
sinom = sin(omega);
idx = 1-cosom>bigeps;       % 2*sin(omega/2)^2>eps, eps1<omega<=pi/2

ti = 2*t.*(1-t);
yy = zeros(ndim, Np);

id = idx;       % eps1<omega<=pi/2
if any(id),
    sclp = sin((1-ti(id)).*omega(id))./sinom(id);
    sclq = sin(ti(id).*omega(id))./sinom(id);
    yy(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

id = ~idx;                % 0<=omega<=eps1
if any(id),
    sclp = 1-ti(id);
    sclq = ti(id);
    yy(:,id) = sclp(nd1,:).*p(:,id)+sclq(nd1,:).*q(:,id);
end;

return;



%% test
ndim = 3;
n = 80;
delta = 1/20;
x1= 1:n;
xx1=1:.01:n;
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
y0 = [zeros(1,n); y1];
yy0 = squad(x1,y0,xx1);
yy1 = yy0(2:4,:);

set(gcf,'color','w');
h = plot3(y1(1,:), y1(2,:), y1(3,:), 'ro', yy1(1,:), yy1(2,:), yy1(3,:), 'b-');
axis equal off;
axis([-1  1 -1 1 -1 1]);

% rotate(h,zdir,90)     % 旋转句柄对象

for i = 1:n,
    text(y1(1,i)+delta, y1(2,i)+delta, y1(3,i)+delta, num2str(i), 'color', 'r');
end;

filename = 'squad.gif';
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
