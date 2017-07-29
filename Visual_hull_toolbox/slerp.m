function yy = slerp(x,y,xx)
% SLERP - Spherical Linear interpolation (within spans less than pi).
%
%   YY = SLERP(X,Y,XX) provide NDIM(>=2) dimensional spherical linear
%   interpolation YY at the values of XX. X are monotonic parameters like
%   time -- N stops. Y(1:NDIM, 1:N) are assumed to be N unit vectors on the
%   NDIM-D sphere. XX are the positions when to interpolate.
%
%   Algorithm: Advanced Animation And Rendering Techniques - Theory And Practice, P364.
%   See also slerp90, spline_interp1, rodrigues, trans_quat_axis, trans_quat_mat.

% by zpf, form BIT, 2015-7-12

bigeps = 1e-12;   % sqrt(bigeps) = 1e-6

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
h = diff(x);
assert(all(h>0) || all(h<0), 'The first argument must be monotonic!');

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

p = y(:,1:n-1);
q = y(:,2:n);
cosom = sum(p.*q);               % 0<=omega<=pi, 0=<omega/2<=pi/2
id1 = 1+cosom>bigeps;       % 2*cos(omega/2)^2 > eps, 0=<omega/2<pi/2-eps1
id2 = 1-cosom>bigeps;       % 2*sin(omega/2)^2 > eps, eps1<omega/2<=pi/2
omega = acos(cosom);
sinom = sin(omega);

assert(isvector(xx),'The third argument must be a vector!');
xx = xx(:)';
Np = length(xx);
[~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).

% now go to local coordinates
t = (xx-x(index))./h(index);
% local values
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
    t = t(id);
    yyy = zeros(ndim, sum(id));
    
    idx = t<tom1;
    ty = t(idx)./tom1(idx);
    sinom = sin(om1);
    sclp = sin((1-ty).*om1(idx))./sinom(idx);
    sclq = sin(ty.*om1(idx))./sinom(idx);
    yyy(:,idx) = sclp(nd1,:).*pl(:,idx) + sclq(nd1,:).*y(:,idx);
    
    idx = ~idx;
    ty = (t(idx)-tom1(idx))./tom2(idx);
    sinom = sin(om2);
    sclp = sin((1-ty).*om2(idx))./sinom(idx);
    sclq = sin(ty.*om2(idx))./sinom(idx);
    yyy(:,idx) = sclp(nd1,:).*y(:,idx) + sclq(nd1,:).*ql(:,idx);
    
    yy(:,id) = yyy;
end;

return;



%% test
ndim = 3;
n = 10;
x1= 1:n*2;
xx1=1:0.1:n*2;
y1=randn(ndim,n);
ny1 = sqrt(sum(y1.*y1));
y1 = y1./ny1(ones(ndim,1),:);
y1 = reshape([y1;-y1],ndim,[]);
yy1 = slerp(x1,y1,xx1);
figure(1);
if ndim ==2,
    plot(y1(1,:),y1(2,:),'ro', yy1(1,:),yy1(2,:),'b-');
elseif ndim == 3,
    plot3(y1(1,:),y1(2,:),y1(3,:),'ro', yy1(1,:),yy1(2,:),yy1(3,:),'b-');
end;
axis equal;

%%
ndim = 3;
n = 10;
delta = 1/20;
x1= 1:n;
xx1=1:.02:n;
y1=randn(ndim,n);
ny1 = sqrt(sum(y1.*y1));
y1 = y1./ny1(ones(ndim,1),:);
yy1 = slerp(x1,y1,xx1);
set(gcf,'color','w');
h = plot3(y1(1,:), y1(2,:), y1(3,:), 'ro', yy1(1,:), yy1(2,:), yy1(3,:), 'b-');
axis equal off;
axis([-1  1 -1 1 -1 1]);

% rotate(h,zdir,90)     % 旋转句柄对象

for i = 1:n,
    text(y1(1,i)+delta, y1(2,i)+delta, y1(3,i)+delta, num2str(i), 'color', 'r');
end;

filename = 'slerp.gif';
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
