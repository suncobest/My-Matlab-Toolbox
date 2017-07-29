function Iout = homography_image( In, ROI, wh, bgcolor)
% homography_image remap image In in ROI to Iout.
% ROI is a quadrangle region surrounded by four counter-clockwise points.
% ROI use MATLAB axis convention: the most left top pixel is [1;1] instead of [0;0]
% Algorithm: compute homography

if nargin<4,
    bgcolor = 1;   % white background;
    if nargin<2,
        Iout = In;
        return;
    end;
end;

assert(isequal(size(ROI),[2,4]),'Unexpected dimension for the 2nd argument!');
% treat ROI as convex hull
xi = [diff([ROI,ROI(:,1)],1,2); zeros(1,4)];
xo = cross(xi(:,4),xi(:,1));
hand = sign(xo(3));
for i =1:3,
    xo = cross(xi(:,i),xi(:,i+1));
    assert(hand==sign(xo(3)), 'The region of interest is not a convex hull!');
end;

if nargin<3 || isempty(wh),
    w = round(norm(ROI(:,4)-ROI(:,1)));
    h = round(norm(ROI(:,2)-ROI(:,1)));
else
    n = length(wh);
    if n ==1,
        w = wh;
        h = wh;
    elseif n==2,
        w = wh(1);
        h = wh(2);
    else
        error('Unexpected input for the 3rd argument!');
    end;
end;

if islogical(In),
    bwflag = 1;
else
    bwflag = 0;
end;

[ny, nx, nc] = size(In);
In = im2double(In);

% check if ROI is out of image.
% extend 1 pixel out of the boundingbox of ROI for interpolation
left = round(min(ROI(1,:)))-1.5;
top = round(min(ROI(2,:)))-1.5;
right = round(max(ROI(1,:)))+1.5;
bottom = round(max(ROI(2,:)))+1.5;

left = min(0.5, left);
top = min(0.5, top);
right = max(nx+0.5, right);
bottom = max(ny+0.5, bottom);
width = right-left;
height = bottom-top;
if width>nx || height>ny,
    if bgcolor,
        Iout = ones(height,width,nc);
    else
        Iout = zeros(height,width,nc);
    end;
    % shift of frame
    left = 0.5-left;
    top = 0.5-top;
    ROI = ROI + repmat([left;top],1,4);
    Iout(top+1:top+ny, left+1:left+nx, :) = In;
    In = Iout;
    nx = width;
    ny = height;
end;

[px,py] = meshgrid(1:w, 1:h);
px = px(:)';
py = py(:)';
% vertices of Iout border in counter-clockwise
xo = [0, 0, w, w; 0, h, h, 0]+0.5;
Homo = compute_homography_lm(ROI, xo);
xo = [px; py; ones(1,w*h)];
xi = Homo*xo;
xi = xi(1:2,:)./(ones(2,1)*xi(3,:));
x0 = floor(xi(1,:));
y0 = floor(xi(2,:));
alpha_x = xi(1,:)-x0;
alpha_y = xi(2,:)-y0;

clt = (1 - alpha_y).*(1 - alpha_x);
crt = (1 - alpha_y).*alpha_x;
clb = alpha_y .* (1 - alpha_x);
crb = alpha_y .* alpha_x;

ind_lt = (x0-1)*ny+y0;
ind_rt = x0*ny+y0;
ind_lb = (x0-1)*ny +y0+1;
ind_rb = x0*ny +y0+1;

Iout = ones(h,w,nc);
Io = Iout(:,:,1);
for i=1:nc,
    Ic = In(:,:,i);
    Io(:) =  clt .* Ic(ind_lt) + crt .* Ic(ind_rt) + clb .* Ic(ind_lb) + crb .* Ic(ind_rb);
    Iout(:,:,i) = Io;
end;

if bwflag,
    Iout = (Iout~=0);
else
    Iout = uint8(255*Iout);
end;

return;



%%  test
nx = 200;
ny = 200;
n = 4;

I=imread('baboon200.jpg');
figure(1);
image(I);
if ismatrix(I),
    colormap(gray(256));
end;
set(1,'color',[1 1 1]);
axis tight
hold on;
flag = 1;
while flag,
%     xv=rand(1,n)*nx*2-100;
%     yv=rand(1,n)*ny*2-100;
    
    [xv,yv] = ginput(4);
    xv = xv';
    yv = yv';
    
    flag = 0;
    vx = [diff([xv,xv(1); yv,yv(1)],1,2); zeros(1,4)];
    vz = cross(vx(:,4),vx(:,1));
    h = sign(vz(3));
    for i =1:3,
        vz = cross(vx(:,i),vx(:,i+1));
        if h~=sign(vz(3)),
            flag = 1;
            disp('Not convex hull, please input again!');
            break;
        end;
    end;
end;
plot([xv, xv(1)], [yv, yv(1)],'b-','linewidth',2.0);
for k=1:4,
    text(xv(k),yv(k),num2str(k),'color','y','fontsize',16);
end;
hold off;

Io = homography_image( I, [xv;yv],[],0);
figure(2);
image(Io);
set(2,'color',[1 1 1]);
axis image;
