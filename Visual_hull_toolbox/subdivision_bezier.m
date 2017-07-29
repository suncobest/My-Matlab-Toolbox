function yy = subdivision_bezier(y, nlevel)
% SUBDIVISION_BEZIER - Recursive subdivision of Bezier spline curves.
%
%   YY = SUBDIVISION_BEZIER(Y, NLEVEL) provide NLEVEL of subdivision of
%   Bezier spline curves (Y). Y is assumed to be the control points of
%   3-order Bezier spline curve. The subroutine XX = HALF_BEZIER(X)
%   calculate 1 time subdivision of X. The main routine call the subroutine
%   NLEVEL times and return the control points of refined splines.
%
%   Algorithm: Advanced Animation And Rendering Techniques - Theory And
%   Practice, P79~81.

% by zpf, form BIT, 2015-8-1

[ndim,ny]=size(y);
if ny==1,
    y=y'; ny=ndim; ndim=1;
end;
assert(ny>=4 && mod(ny, 3)==1, 'Unexpected number of control points!');

if nargin<2,
    nlevel = 1;
end;
n = (ny-1)/3;
yn = y(:,ny);
yy = reshape(y(:,1:ny-1), ndim, 3, n);
yy = [yy, cat(3,yy(:,1, 2:n), yn)];

for i=1:nlevel,
    yy = half_bezier(yy, ndim);
end;
yy = [reshape(yy(:,1:3, :), ndim, []), yn];

return;


function xx = half_bezier(x, ndim)
% HALF_BEZIER - Divide the 3-order bezier curve at the half.
q1 = x(:, 1, :);
r4 = x(:, 4, :);
x= (x(:,1:3,:)+x(:,2:4,:))/2;
q2 = x(:, 1, :);
r3 = x(:, 3, :);
x= (x( :,1:2,:)+x(:,2:3,:))/2;
q3 = x(:, 1, :);
r2 = x(:, 2, :);
q4 = (x( :,1,:)+x(:,2,:))/2;
xx = [q1, q2, q3, q4, q4, r2, r3, r4];
xx = reshape(xx,ndim,4,[]);

return;




%% test
nd = 3;
N = 4;
delta = 1/40;
y1=randn(nd,N);
y2 = y1;   % y2=subdivision_bezier(y1,2);
rgb = 'myr';
if nd ==2,
    hold off;
    plot(y1(1,:), y1(2,:), 'bo-'); hold on;
    for i = 1:N,
        if i>1,
            y2=subdivision_bezier(y2);
            plot(y2(1,:), y2(2,:),'o:', 'color', rgb(i-1));
        end;
        text(y1(1,i)+delta, y1(2,i)+delta,  num2str(i), 'color', 'r');
    end;
elseif nd ==3,
    hold off;
    plot3(y1(1,:), y1(2,:), y1(3,:),'bo-'); hold on;
    for i = 1:N,
        if i>1,
            y2=subdivision_bezier(y2);
            plot3(y2(1,:), y2(2,:),y2(3,:),'o:', 'color', rgb(i-1));
        end;
        text(y1(1,i)+delta, y1(2,i)+delta, y1(3,i)+delta,  num2str(i), 'color', 'r');
    end;
end;
set(gcf,'color','w');
axis equal off;
% h = plot(y1(1,:), y1(2,:), 'ro', y1(1,:), y1(2,:), 'g-', y2(1,:), y2(2,:),'mo', y2(1,:),y2(2,:), 'b-');
% rotate(h,zdir,90)     % Ðý×ª¾ä±ú¶ÔÏó

