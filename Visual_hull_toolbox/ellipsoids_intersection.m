function [intp, ind] = ellipsoids_intersection( cc, ax, vec, XX )
% ELLIPSOIDS_INTERSECTION compute the volume of two ellipsoids intersection in fraction of each
% ellipsoid.
%   cc: (input) center of two ellipsoids. (M-dimension)
%   ax: (input) semi axes of ellipsoid.
%   vec: (input) orientation of ellipsoid axes.
%   XX: (input) points of dimension m*npts; (m>=1), assumed as 1D points if XX is vector;
%   intp: (output) intersection over whole volume inside every ellipsoid. (1*2)
%   ind: (output) index of points inside of data XX in every ellipsoid. (npts*2)
%  See also principal_ellipsoid, optimal_ellipsoid, ellipsoid_fitting.

[m,n] = size(cc);
assert(n==2,'There must be two ellipsoids!');
[m1,n] = size(ax);
assert(n==2 && m1==m,'Unexpected dimension of the 2nd input!');
[m1,n] = size(vec);
assert(m1==m && n==m*2,'Unexpected dimension of the 3rd input!');

if nargin<4 || isempty(XX),
    npts = 1e6;
    XX = 2*rand(m,npts)-ones(m,npts);
    XX = XX(:, sum(XX.^2,1)<=1);
    npts = size(XX,2);
    XX = vec(:,1:m).*(ones(m,1)*ax(:,1)')*XX+cc(:,1)*ones(1,npts);
    YY = 1./ax(:,2)*ones(1,m).*vec(:,m+1:end)'*(XX-cc(:,2)*ones(1,npts));
    nn = sum(sum(YY.^2,1)<=1);
    intp = nn/npts*[1, prod(ax(:,2))/prod(ax(:,1))];
    ind = [];
    return;
end;

[m1,npts] = size(XX);
if npts==1,
    npts = m1;
    m1 = 1;
end;
assert(m1==m && npts>0,'Unexpected dimension of the 4th input!');
ind = false(npts,2);
YY = 1./ax(:,1)*ones(1,m).*vec(:,1:m)'*(XX-cc(:,1)*ones(1,npts));
ind(:,1) = sum(YY.^2,1)<=1;
YY = 1./ax(:,2)*ones(1,m).*vec(:,m+1:end)'*(XX-cc(:,2)*ones(1,npts));
ind(:,2) = sum(YY.^2,1)<=1;
intp = sum(ind(:,1) & ind(:,2))./sum(ind,1);

return;



%% Test
% 1D
c = randi(10,1,2);
a= randi(4,1,2);
v = [1 1];
x = -1:0.1:10;
[p,id] = ellipsoids_intersection(c, a, v, x);

%% 2D
mk = 30;
c = rand(2,2)*10;
a= randi(5,2,2);
alpha = randn(1,2);
v =[cos(alpha(1)),-sin(alpha(1)),cos(alpha(2)),-sin(alpha(2)); sin(alpha(1)),cos(alpha(1)), ...
    sin(alpha(2)),cos(alpha(2))] ;
x = rand(2,1000)*10;
[p,id] = ellipsoids_intersection(c, a, v, x);
figure;
hold on;
id1 = id(:,1) & ~id(:,2);
id2 = id(:,2) & ~id(:,1);
id3 = id(:,1) & id(:,2);
id4 = ~(id(:,1) | id(:,2));
rgb = 'rgbc';
for i=1:4,
    eval(['plot(x(1,id' num2str(i) '), x(2,id' num2str(i) '), ''' rgb(i) '.'')']);
end;
t = linspace(0,2*pi,mk);
for i=1:2,
    xy = v(:,(i-1)*2+1:(i-1)*2+2).*(ones(2,1)*a(:,i)')*[cos(t); sin(t)]+c(:,i)*ones(1,mk);
    plot(xy(1,:),xy(2,:),'b-');
end;
axis equal tight off;

%% 3D
mk = 30;
c = rand(3,2)*5;
a= randi(3,3,2);
om = randn(3,2);
v = [rodrigues(om(:,1)),rodrigues(om(:,2))];
x = rand(3,10000)*5;
[p,id] = ellipsoids_intersection(c, a, v, x);
figure;
hold on;
id1 = id(:,1) & ~id(:,2);
id2 = id(:,2) & ~id(:,1);
id3 = id(:,1) & id(:,2);
id4 = ~(id(:,1) | id(:,2));
rgb = 'rgbc';
for i=1:4,
    eval(['plot3(x(1,id' num2str(i) '), x(2,id' num2str(i) '), x(3,id' num2str(i) '), ''' rgb(i) '.'')']);
end;
for i=1:2,
    [X,Y,Z]=ellipsoid(0,0,0,a(1,i),a(2,i),a(3,i),mk);
    xyz = v(:,(i-1)*3+1:(i-1)*3+3)*[X(:)';Y(:)';Z(:)']+c(:,i)*ones(1,(mk+1)^2);
    X(:) = xyz(1,:);
    Y(:) = xyz(2,:);
    Z(:) = xyz(3,:);
    surf(X,Y,Z,'FaceColor', 'b','EdgeColor','g','FaceAlpha',0.01);
end;
axis equal tight off;

