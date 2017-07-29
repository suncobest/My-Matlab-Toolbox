function yy = spline_bezier(x,y,xx)
% SPLINE_BEZIER - Bezier spline interpolation with C1 continuity.
%   YY = SPLINE_BEZIER(X,Y,XX) provide interpolation YY at the values of XX. 
%   See also spline_interp2, spline_interp3.

% by zpf, form BIT, 2015-7-29

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
assert(isvector(xx),'The third argument must be a vector!');
m = size(xx,2);
if m==1,
    xx = xx';
end;

% histogram filter, if [n,bin] = histc(x,edges), n(k) = sum(bin==k).
[~,index] = histc(xx,[-inf,x(2:n-1),inf]);    % n slots with the last slot 0 elements (inf)
nd1 = ones(ndim,1);

% caculate the auxiliary points q
q = y;
t = h(2:n-1)./(2*(h(1:n-2)+h(2:n-1)));
yy = 0.5-t;
q(:, 2:n-1) = y(:,2:n-1)*3/2-t(nd1,:).*y(:,1:n-2)-yy(nd1,:).*y(:,3:n);

% now go to local coordinates
t = (xx-x(index))./h(index);
p1 = y(:,1:n-1);
p2 = y(:,2:n);
q1 = q(:,1:n-1);
q2 = q(:,2:n);
p1 = p1(:, index);
p2 = p2(:, index);
q1 = q1(:, index);
q2 = q2(:, index);

q = 1-t;
p1 = q(nd1,:).*p1+t(nd1,:).*p2;
p2 = q(nd1,:).*q1+t(nd1,:).*q2;
t = 2*t.*q;
q = 1-t;
yy = q(nd1,:).*p1+t(nd1,:).*p2;
if  ndim==1 && m==1,
    yy = yy';
end;


return;



%% test
Np = 5;
Nd = 3;
x1=1:Np;
n2 = 200;
lw = 1.5;
xx1=linspace(1, Np, n2);
y1 =10*rand(Nd,Np);
yy1 = spline_bezier(x1,y1,xx1);

figure(1);
if Nd==2,
    plot(y1(1,:),y1(2,:),'bo',  yy1(1,:),yy1(2,:),'b-','linewidth', lw);
elseif Nd==3,
    plot3(y1(1,:),y1(2,:),y1(3,:), 'bo',  yy1(1,:),yy1(2,:),yy1(3,:),'b-', 'linewidth', lw);
end;
hold on;

kinds = 'NCKQ';
rgb = 'gckr';
for i=1:4,
    yy1 = spline_interp3(x1,y1,xx1,kinds(i));
    figure(1);
    if Nd==2,
        plot(yy1(1,:),yy1(2,:),'color', rgb(i), 'linewidth', lw);
    elseif Nd==3,
        plot3(yy1(1,:),yy1(2,:),yy1(3,:),'color', rgb(i), 'linewidth', lw);
    end;
end;

yy1 = spline_interp2(x1,[zeros(Nd,1), y1],xx1);
figure(1);
if Nd==2,
    plot(yy1(1,:),yy1(2,:),'m-', 'linewidth', lw);
elseif Nd==3,
    plot3(yy1(1,:),yy1(2,:),yy1(3,:),'m-', 'linewidth', lw);
end;

yy1 = spline_interp2(x1, y1, xx1,1);
figure(1);
if Nd==2,
    plot(yy1(1,:),yy1(2,:),'y-', 'linewidth', lw);
elseif Nd==3,
    plot3(yy1(1,:),yy1(2,:),yy1(3,:),'y-', 'linewidth', lw);
end;

axis equal;

%% see draw_spline_bezier
