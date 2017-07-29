% calibration
if ~exist('imname','var') || ~exist('I','var')
    read_image;
end

figure(2);
image(I);
colormap(gray(256));
set(2,'color',[1 1 1]);
axis image;

% set the size of window for finding corners
wint_default = 7;
disp('Window size for corner finder (wintx and winty):');
wintx = input(['wintx ([] = ' num2str(wint_default) ') = ']);
if isempty(wintx), wintx = wint_default; end;
wintx = round(wintx);
winty = input(['winty ([] = ' num2str(wint_default) ') = ']);
if isempty(winty), winty = wint_default; end;
winty = round(winty);
fprintf(1,'Window size = %d x %d\n',2*wintx+1,2*winty+1);


title('Click on the four vertexes of a known rectangle (first vertex = origin)');
disp('Click on the four vertexes of a known rectangle (first vertex = origin)...');

x= [];y = [];
figure(2); hold on;
for count = 1:4
    [xi,yi] = ginput(1);
    xxi = cornerfinder3([xi;yi],I,winty,wintx);
    xi = xxi(1);
    yi = xxi(2);
    x = [x;xi];
    y = [y;yi];
    figure(2);
    plot(xi,yi,'g+','linewidth',2,'Markersize',10);   
    plot(xi + [wintx+.5, -(wintx+.5), -(wintx+.5), wintx+.5, wintx+.5],yi + [winty+.5, winty+.5, -(winty+.5), -(winty+.5), winty+.5],'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);  
    drawnow;
end
drawnow;
hold off;

oxy = cornerfinder3([x';y'],I,winty,wintx); % the four corners

x = oxy(1,:)';
y = oxy(2,:)';

% Sort the corners:
x_mean = mean(x);
y_mean = mean(y);
x_v = x - x_mean;
y_v = y - y_mean;

theta = atan2(-y_v,x_v);  % atan(-y,x)=-atan(y,x)； 

% 若theta = atan2(y_v,x_v)对应的幅角为(a1,a2,a3,a4), 得到角度差theta-theta(1)为(a1-a1,a2-a1,a3-a1,a4-a1)=angle
% 则theta = atan2(-y_v,x_v)的幅角度为(-a1,-a2,-a3,-a4)，角度差theta-theta(1)为(a1-a1,a1-a2,a1-a3,a1-a4)=-angle
% 在图像坐标系中，幅角的正方向(从x轴旋转到矢量)为顺时针angle，法线朝里；但棋盘图案的世界坐标系的xy平面在图像上的投影逆时针为正(-angle)，法线朝外。
% mod(x,y)=x-n.*y 其中 n=floor(x./y)。由于角度差a满足-2*pi<=a<=2*pi，若a>=0,则mod(a,2*pi)=a；若a<0,则mod(a,2*pi)=2*pi+a。
% 于是0<=mod(a,2*pi)<2*pi，将a调整到[0,2*pi)。对-angle进行排序即是从原点逆时针排序。

[junk,ind] = sort(mod(theta-theta(1),2*pi));  % ind: counterclockwise ascending order; first point is the orgin (o,x,xy,y)

x = x(ind);
y = y(ind);
x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);   % (x1,y1) is the origin, (x2,y2) is the extreme point in x direction
y1= y(1); y2 = y(2); y3 = y(3); y4 = y(4);   % (x3,y3) is the extreme point in the diagonal direction, (x4,y4) is the extreme point in y direction

oxy = oxy(:,ind);

% Find center:
p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
x5 = p_center(1)/p_center(3);
y5 = p_center(2)/p_center(3);

% center on the X axis:
x6 = (x1 + x2)/2;
y6 = (y1 + y2)/2;

% center on the Y axis:
x7 = (x1 + x4)/2;
y7 = (y1 + y4)/2;

% Direction of displacement for the X axis:
vX = [x6-x5;y6-y5];
vX = vX / norm(vX);

% Direction of displacement for the Y axis:
vY = [x7-x5;y7-y5];
vY = vY / norm(vY);

% Direction of displacement for the origin:
vO = [x1 - x5; y1 - y5];
vO = vO / norm(vO);

delta = min(nx,ny)/30;

figure(2); 
image(I);
colormap(gray(256));
set(2,'color',[1 1 1]);
axis image;
hold on;
plot(x,y,'g+','linewidth',2,'Markersize',10);
plot([x;x(1)],[y;y(1)],'color',[ 1.00 0.50 0.30 ],'linewidth',2);
arrow3(oxy(:,1)',oxy(:,2)','b2',1.2,2);
arrow3(oxy(:,1)',oxy(:,4)','b2',1.2,2);
hx=text('Interpreter','latex','String','$$x$$','Position',[x6 + delta * vX(1), y6 + delta*vX(2)]);
set(hx,'color','g','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
hy=text('Interpreter','latex','String','$$y$$','Position',[x7 + delta * vY(1), y7 + delta*vY(2)]);
set(hy,'color','g','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
hO=text('Interpreter','latex','String','$$o$$','Position',[x1 + delta * vO(1) ,y1 + delta*vO(2)]);
set(hO,'color','g','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

hold off;


% set the parameters for calibration

% focal length
disp('Focal length of camera:');
focus = input('focus ([] = 1) = ');
if isempty(focus), focus = 1; end;
fprintf(1,'Focal length:    focus = %3.5f\n',focus);


% pixel size
disp('Pixel size of the camera CCD:');
pix = input('pix ([] = 1) = ');
if isempty(pix), pix = 1; end;
fprintf(1,'pixel size:    pix = %3.5f\n',pix);


% physical length of the two edges
disp('Physical length of the two perpendicular edges (lx and ly):');
lx = input('length of lx ([] = 1) = ');
if isempty(lx), lx = 1; end;
ly = input('length of ly ([] = 1) = ');
if isempty(ly), ly = 1; end;
fprintf(1,'Length of lx:   lx = %3.5f \nLength of ly:   ly= %3.5f\n',lx,ly);


% the sign of Z component
disp('Signs of Z components of the two edges (z_sign_lx and z_sign_ly): please input 1 or -1');
z_sign_lx = input('Z sign of lx ([] = 1) = ');
if isempty(z_sign_lx), z_sign_lx = 1; end;
z_sign_ly = input('Z sign of ly ([] = 1) = ');
if isempty(z_sign_ly), z_sign_ly = 1; end;
fprintf(1,'z_sign_lx = %d \nz_sign_ly = %d\n',z_sign_lx,z_sign_ly);

alpha = 90; 

% Calculate the length of the two edges (in pixels) and the angle between them
pedge1 = [x2-x1, y2-y1];
pedge2 = [x4-x1, y4-y1];
plx = norm(pedge1);
ply = norm(pedge2);
cbeta = sum(pedge1.*pedge2)/(plx*ply);
beta = acosd(cbeta);

rx = plx/lx;
ry = ply/ly;

% Calibration 
% magnification
mag = pix/sind(alpha) * sqrt( ( rx^2 + ry^2 - 2*cosd(alpha)*cbeta*rx*ry + ...
      sqrt( rx^4 + ry^4 + 2*(cosd(2*alpha) + cosd(2*beta) +1)*rx^2*ry^2 - ...
      4*cosd(alpha)*cbeta*rx*ry*(rx^2 + ry^2) ) ) /2);

% distance from object to optical center: so
% distance from image to optical center: si
so = focus/mag + focus;
si = so*mag;

% Vectors of the two edges in the scene
edge1 = [pedge1*pix/mag, z_sign_lx * sqrt(lx^2-(plx*pix/mag)^2)];
edge2 = [pedge2*pix/mag, z_sign_ly * sqrt(ly^2-(ply*pix/mag)^2)];

% direction of the plane in the scene
dir_plane = cross(edge1,edge2);
dir_plane = dir_plane/norm(dir_plane);
theta = acosd(dir_plane(3));

% compute the homography from the plane of interest to the image
% 由于[x;y;1] = Homo*[0 1 1 0;0 0 1 1;1 1 1 1]，且[0 1 1 0;0 0 1 1;1 1 1 1]=diag([1/lx,1/ly,1])*[0 lx lx 0;0 0 ly ly;1 1 1 1]
% 所以[x;y;1] = Homo*diag([1/lx,1/ly,1])*[0 lx lx 0;0 0 ly ly;1 1 1 1]
Homo = compute_homography_lm([x y ones(4,1)]',[0 1 1 0;0 0 1 1;1 1 1 1]);  % (o,x,xy,y)
Homo = Homo*diag([1/lx,1/ly,1]);


disp('Plane calibrated!');
fprintf(1,'\nCalibration result:\n--------------------------------------------------------------------------------\n');
fprintf(1,'Magnification of the image:                   mag = %3.5f\n',mag);
fprintf(1,'Distance from object to the optical center:    so = %3.5f\n',so);
fprintf(1,'Distance from image to the optical center:     si = %3.5f\n',si);
fprintf(1,'Direction of the calibrated plane:      dir_plane = [%3.5f, %3.5f, %3.5f]\n',dir_plane);
fprintf(1,'Angle between the plane and the image:      theta = %3.5f in degree\n',theta);
fprintf(1,'--------------------------------------------------------------------------------\n');

