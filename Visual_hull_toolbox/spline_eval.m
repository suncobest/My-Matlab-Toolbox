function [yy, dyydxx, d2yyd2xx] = spline_eval(x, coefs, xx)
% SPLINE_EVAL:  Evaluate piecewise polynomial of order 1~3.
%   YY = spline_eval(X, COEFS, XX) returns the value, at the entries of XX, of the 
%   piecewise polynomial coefficients in COEFS, as constructed by spline_interp1, 
%   spline_interp2 and spline_interp1.
%   COEFS((j-1)*ndim+1:j*ndim,:) is the j-th polynomial of the spline curve.
%   If size(coefs,2)==4:
%   coefficient [a_i,b_i,c_i,d_i]:   s_i = a_i*t_i^3+b_i*t_i^2+c_i*t_i+d_i, t_i = (x-x_i)/h(i);

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
n = length(x);
assert(n>=2,'It takes at least two points to make a spline curve!');
h = diff(x);
assert(all(h>0) || all(h<0),'The first argument must be monotonic!');
[mc, nc] = size(coefs);
ndim = mc/(n-1);
assert(ndim==floor(ndim),'Unexpected rows of the 2nd argument!');
if h(1)<0,
    n1 = n:-1:1;
    x = x(n1);
    h = -h(n1(2:end));
end;

nd1 = ones(ndim,1);
if nargin<3,
    xx= x;
    Np = n;
else
    assert(isvector(xx),'The third argument must be a vector!');
    xx = xx(:)';
    Np = length(xx);
end;
[~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).
% now go to local coordinates
t = (xx-x(index))./h(index);

if ndim>1,  % replicate t and index in case Points are vector
    t = reshape(t(nd1,:),Np*ndim,1);
    index = ndim*index; temp = (-ndim:-1)';
    index = reshape(1+index(nd1,:)+temp(:,ones(1,Np)), Np*ndim,1);
end;
v = coefs(index,1);
for i = 2:nc,
    v = v.*t(:)+coefs(index,i);
end;
yy = reshape(v,ndim,Np);

if nargout>1,
    nm1 = ones(mc,1);
    temp = nc-1:-1:1;
    coefs = coefs(:,1:nc-1).*temp(nm1,:)./(reshape(h(nd1,:),mc,1)*ones(1,nc-1));
    v = coefs(index,1);
    for i = 2:nc-1,
        v = v.*t(:)+coefs(index,i);
    end;
    dyydxx = reshape(v,ndim,Np);
    if nargout>2,
        temp = temp(2:end);
        coefs = coefs(:,1:nc-2).*temp(nm1,:)./(reshape(h(nd1,:),mc,1)*ones(1,nc-2));
        v = coefs(index,1);
        for i = 2:nc-2,
            v = v.*t(:)+coefs(index,i);
        end;
        d2yyd2xx = reshape(v,ndim,Np);
    end;
end;

return;


%% Test
theta = [0 30 45 60 90 120 135 150 180]/180*pi;
r = [cos(theta); sin(theta)];
th=(0:5:180)/180*pi;
[r1, cc]= spline_interp3(theta,r,th);
plot(r(1,:),r(2,:),'ro', r1(1,:),r1(2,:),'g.');
axis equal;
[r2, dr2, d2r2] = spline_eval(theta, cc, th);
err = norm(r1-r2)
err = norm(r1-[cos(th); sin(th)])
err = norm([-sin(th); cos(th)] - dr2)
err = norm([-cos(th); -sin(th)] - d2r2)
[r1, dr1, d2r1] = spline_eval(theta, cc);
err = norm(r-r1)
dr=[-sin(theta); cos(theta)], dr1
d2r=[-cos(theta); -sin(theta)], d2r1

%%
x = sort([-rand(1,3),rand(1,3)]);
y = [x.^2; x.^3; exp(x)];
xx = -1:0.05:1;
[yy, cc] = spline_interp3(x,y,xx);
plot3(y(1,:),y(2,:),y(3,:),'ro', yy(1,:),yy(2,:),yy(3,:),'g.');
axis equal;
[yy1, dyy1, d2yy1] = spline_eval(x, cc, xx);
err = norm(yy-yy1)
err = norm([xx.^2; xx.^3; exp(xx)]-yy1)
err = norm([2*xx; 3*xx.^2; exp(xx)]-dyy1)
err = norm([2*ones(1,size(xx,2)); 6*xx; exp(xx)]-d2yy1)
[y1, dy1, d2y1] = spline_eval(x, cc);
dy=[2*x; 3*x.^2; exp(x)], dy1
d2y=[2*ones(1,size(x,2)); 6*x; exp(x)], d2y1
