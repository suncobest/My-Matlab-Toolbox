function [yy,coefs] = spline_interp1(x,y,xx)
% SPLINE_INTERP1 linear spline interpolation
%   COEFS = SPLINE_INTERP1(X,Y) provides the piecewise polynomial form of the 
%   linear spline curve to interpolate the values Y at the break sites X.
%   If Y is a vector, then Y(j) is taken as the value to be matched at X(j), 
%   hence Y must be of the same length as X.
%   If Y is a matrix, then Y(:,j) is taken as the value to be matched at X(j),
%   hence the columns of Y must equal length(X).
%   COEFS((j-1)*ndim+1:j*ndim,:) is the j-th polynomial of the spline curve.
%   coefficient [a_i,b_i]: s_i = a_i*t_i+b_i, ti = (x-x_i)/h(i);
%
%   YY = SPLINE_INTERP1(X,Y,XX) provide interpolation YY at the values of XX. 
%   See also spline_interp2, spline_interp3.

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
assert(n==ny,'Number of breaks do not match!');
if h(1)<0,
    n1 = n:-1:1;
    x = x(n1);
    y = y(:,n1);
    h = -h(n1(2:end));
end;

b = y(:,1:n-1);
a = y(:,2:n)-b;

coefs =[a(:),b(:)];
if nargin==2,
    yy = coefs;
elseif nargin==3,
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
        nd1 = ones(ndim,1);
        t = reshape(t(nd1,:),Np*ndim,1);
        index = ndim*index; temp = (-ndim:-1)';
        index = reshape(1+index(nd1,:)+temp(:,ones(1,Np)), Np*ndim,1);
    end;
       
    v = coefs(index,1).*t(:)+coefs(index,2); 
    
    yy = reshape(v,ndim,Np);   
    if  ndim==1 && m==1,
        yy = yy';      
    end;
end;

return;



%% test

x1=-1:.2:1;
xx1=-2:.1:2;
y1=x1.^4;
z1=x1.^3;
yy1 = spline_interp1(x1,[y1;z1],xx1);
plot(x1,y1,'ro', x1,z1,'yo', xx1,yy1(1,:),'b-', xx1,yy1(2,:),'g-');
