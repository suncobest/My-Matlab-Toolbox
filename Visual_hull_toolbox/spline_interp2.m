function [yy,coefs] = spline_interp2(x,y,xx, SW)
% SPLINE_INTERP2 - quadratic spline interpolation
%   COEFS = SPLINE_INTERP2(X,Y) provides the piecewise polynomial form of the
%   quadratic spline curve to interpolate the values Y at the break sites X.
%   If Y is a vector, then Y(j) is taken as the value to be matched at X(j),
%   hence Y must be of the same length as X.
%   If Y is a matrix, then Y(:,j) is taken as the value to be matched at X(j),
%   hence the columns of Y must equal length(X).
%   COEFS((j-1)*ndim+1:j*ndim,:) is the j-th polynomial of the spline curve.
%   coefficient [a_i,b_i,c_i]: s_i = a_i*t_i^2+b_i*t_i+c_i, t_i = (x-x_i)/h(i);
%   See the exception in case 2.
%   Let n = length(x); if n==2, then the curve is degenerated into linear spline;
%
%   SW is the switch to turn on case 1 (and turn off case 3) if SW==1 (default=0), .
%
%   Three cases: (Ordinarily, the not-a-knot end conditions (case 3) are used.)
%   Case1, the first segment is a straight line, so a1=0, b1=y2-y1; in this case,
%   Y must be of the same length as X;
%
%   Case 2, given the tangent of the 1st knot, so s1_x(0)=b1/h1=B1; in this case,
%   Y must add a tangent element before the 1st knot.
%
%   Case 3, the 1st and 2nd pieces of curve is the same quadratic curve. In other
%   words, the 2nd knot is not a knot any more. b_1/h2(1)= b_2/h2(2), 
%   Y must be of the same length as X in this case;
%
%   YY = SPLINE_INTERP2(X,Y,XX) provide quatratic spline interpolation YY at the values of XX.
%
%   See also spline_interp1, spline_interp3.

% by zpf, form BIT, 2015-7-1

assert(isvector(x),'The first argument must be a vector!');
x = x(:)';
n = length(x);
assert(n>=2,'It takes at least two points to make a spline curve!');
h = diff(x);
assert(all(h>0) || all(h<0),'The first argument must be monotonic!');

[ndim,ny]=size(y);
if ny==1,
    y=y'; ny=ndim; ndim=1;
end;
if h(1)<0,
    n1 = n:-1:1;
    x = x(n1);
    h = -h(n1(2:end));
    y = y(:,ny:-1:1);
end;

if nargin<4,
    SW = false;
else
    SW = ~~SW;
end;
b = zeros(ndim,n-1);
nd1 = ones(ndim,1);

if n==ny,
    dy = diff(y,[],2);
    delta = dy ./ h(nd1,:);
    c = y(:,1:n-1);
    if n==2 || SW,   % Case 1
        b(:,1) = dy(:,1);
    else                    % Case 3
        a = (delta(:,1)*h(2)^2+dy(:,2)*h(1))/(h(1)+h(2));
        b(:,1) = 2*dy(:,1)-a*h(1)/h(2);
    end;
elseif ny==n+1      % Case 2
    dy = diff(y(:,2:n+1),[],2);
    delta = dy ./ h(nd1,:);
    b(:,1) = y(:,1)*h(1);
    c = y(:,2:n);
else
    error('Number of breaks do not match!');
end;

for i = 2:n-1,
    b(:,i) = (2*delta(:,i-1)-b(:,i-1)/h(i-1))*h(i);
end;
a = dy - b;

% coefficients of splines
coefs =[a(:),b(:),c(:)];

if nargin==2,
    yy = coefs;
elseif nargin>=3,
    assert(isvector(xx),'The third argument must be a vector!');
    m = size(xx,2);
    if m==1,
        xx = xx';
    end;
    Np = length(xx);
    [~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).
    % now go to local coordinates
    t = (xx-x(index))./h(index);
    
    if ndim>1,  % replicate t and index in case Points are vector
        t = reshape(t(nd1,:),Np*ndim,1);
        index = ndim*index; temp = (-ndim:-1)';
        index = reshape(1+index(nd1,:)+temp(:,ones(1,Np)), Np*ndim,1);
    end;
    
    v = coefs(index,1).*t(:)+coefs(index,2);
    v = v.*t(:)+coefs(index,3);
    
    yy = reshape(v,ndim,Np);
    if  ndim==1 && m==1,
        yy = yy';
    end;
end;

return;



%% test
n1 = 8;
n2 = 400;
interval  = [-12 12 -5  5];
x1=linspace(-10, 10, n1);
xx1=linspace(-11, 11, n2);
y1=sin(x1);
z1=cos(x1.^3);

[yy1,coef] = spline_interp2(x1,[y1;z1],xx1);
figure(1);
subplot(2,2,1);
plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
title('Not-a-knot start');
axis(interval);

[yy1,coef]  = spline_interp2(x1, [y1;z1], xx1,1);
subplot(2,2,2);
plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
title('Straight-line start');
axis(interval);

[yy1,coef]  = spline_interp2(x1, [0,y1;0,z1],xx1);
subplot(2,2,3);
plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
title('Zero-tangent start');
axis(interval);

[yy1,coef]  = spline_interp2(x1, [1,y1;-0.5,z1],xx1);
subplot(2,2,4);
plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
title('Non-zero-tangent start');
axis(interval);



