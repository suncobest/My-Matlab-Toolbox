function IN = region_inpolygon(nx, ny, xv,yv)
% This function will output a binary image containing a polygon region
% given the size of image (nx, ny), and corner points.
% Image axis [XMIN XMAX YMIN YMAX] = [0, nx, 0, ny]+0.5;
% The center of pixel of i row j col is at the position (j, i).
% See also inpolygon, convex_hull, convhull, convhulln.

% generate meshgrid for inpolygon
[X,Y] = meshgrid(1:nx, 1:ny);
IN = inpolygon(X,Y, xv, yv);
return;


%% Test
nx = 400;
ny = 300;
n = 6;
xv=rand(1,n)*nx;
yv=rand(1,n)*ny;
tic;
I=region_inpolygon(nx,ny,xv,yv);
t=toc
figure(1);
image(I);
colormap(gray(2));
set(1,'color',[1 1 1]);
axis image;
hold on;
plot([xv, xv(1)], [yv, yv(1)],'go',[xv, xv(1)], [yv, yv(1)],'b-','linewidth',2.0);
hold off;