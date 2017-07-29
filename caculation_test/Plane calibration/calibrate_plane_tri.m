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


title('Click on the three vertexes of a known angle in the counter-clockwise direction');
disp('Click on the three vertexes of a known angle in the counter-clockwise direction (the second vertex is the corner)...');

x= [];y = [];
figure(2); hold on;
for count = 1:3
    [xi,yi] = ginput(1);
    xxi = cornerfinder3([xi;yi],I,winty,wintx);
    xi = xxi(1);
    yi = xxi(2);
    x = [x;xi];
    y = [y;yi];
    figure(2);
    if count>1
        plot(x(count-1:count),y(count-1:count),'b-','linewidth',2.5);
    end
    plot(xi,yi,'g+','linewidth',2,'Markersize',10);   
    plot(xi + [wintx+.5, -(wintx+.5), -(wintx+.5), wintx+.5, wintx+.5],yi + [winty+.5, winty+.5, -(winty+.5), -(winty+.5), winty+.5],'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);  
    drawnow;
end
hold off;

% Calculate the length of the two edges (in pixels) and the angle between them
pedge1 = [x(1)-x(2), y(1)-y(2)];
pedge2 = [x(3)-x(2), y(3)-y(2)];
d1 = norm(pedge1);
d2 = norm(pedge2);
cbeta = sum(pedge1.*pedge2)/(d1*d2);
beta = acosd(cbeta);


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


delta = min(nx,ny)/30;
radius = 3*delta;
dstart = atan2d(pedge1(2),pedge1(1));
dend = atan2d(pedge2(2),pedge2(1));

if dend-dstart>180,
    dend=dend-360;
end;

position_angle = [x(2),y(2)] + (delta+radius)*[cosd((dstart+dend)/2), sind((dstart+dend)/2)];

% draw the two vectors and make annotations
AX=figure(2);
image(I);
colormap(gray(256));
set(AX,'color',[1 1 1]);
axis image;
title('Set the parameters for calibration');
hold on;
arrow3([x(2),y(2)],[x(1),y(1)],'b2',1.2,2);
arrow3([x(2),y(2)],[x(3),y(3)],'b2',1.2,2);
draw_arc(AX,[x(2),y(2)],radius,dstart,dend,20);
position_l1 = [x(1),y(1)] + delta*pedge1/d1;
position_l2 = [x(3),y(3)] + delta*pedge2/d2;
text('Interpreter','latex','String','$$l_1$$','Position',position_l1,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
text('Interpreter','latex','String','$$l_2$$','Position',position_l2,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
text('Interpreter','latex','String','$$\alpha$$','Position',position_angle,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
hold off;

% physical length of the two edges
disp('Physical length of the two edges (l1 and l2):');
l1 = input('length of l1 ([] = 1) = ');
if isempty(l1), l1 = 1; end;
l2 = input('length of l2 ([] = 1) = ');
if isempty(l2), l2 = 1; end;
fprintf(1,'Length of l1:   l1 = %3.5f \nLength of l2:   l2= %3.5f\n',l1,l2);


% the sign of Z component
disp('Signs of Z components of the two edges (z_sign_l1 and z_sign_l2): please input 1 or -1');
z_sign_l1 = input('Z sign of l1 ([] = 1) = ');
if isempty(z_sign_l1), z_sign_l1 = 1; end;
z_sign_l2 = input('Z sign of l2 ([] = 1) = ');
if isempty(z_sign_l2), z_sign_l2 = 1; end;
fprintf(1,'z_sign_l1 = %d \nz_sign_l2 = %d\n',z_sign_l1,z_sign_l2);


% physical angle between the two edges
disp('Physical angle between the two edges (in degree):');
alpha = input('Physical angle alpha ([] = 90) = ');
if isempty(alpha), alpha = 90; end;
fprintf(1,'Physical angle:    alpha = %3.5f\n',alpha);


dl1 = d1/l1;
dl2 = d2/l2;

% Calibration 
% magnification
mag = pix/sind(alpha) * sqrt( ( dl1^2 + dl2^2 - 2*cosd(alpha)*cbeta*dl1*dl2 + ...
      sqrt( dl1^4 + dl2^4 + 2*(cosd(2*alpha) + cosd(2*beta) +1)*dl1^2*dl2^2 - ...
      4*cosd(alpha)*cbeta*dl1*dl2*(dl1^2 + dl2^2) ) ) /2);

% distance from object to optical center: so
% distance from image to optical center: si
so = focus/mag + focus;
si = so*mag;

% Vectors of the two edges in the scene
edge1 = [pedge1*pix/mag,z_sign_l1 * sqrt(l1^2-(d1*pix/mag)^2)];
edge2 = [pedge2*pix/mag,z_sign_l2 * sqrt(l2^2-(d2*pix/mag)^2)];

% direction of the plane in the scene
dir_plane = cross(edge1,edge2);
dir_plane = dir_plane/norm(dir_plane);
theta = acosd(dir_plane(3));

disp('Plane calibrated!');
fprintf(1,'\nCalibration result:\n--------------------------------------------------------------------------------\n');
fprintf(1,'Magnification of the image:                   mag = %3.5f\n',mag);
fprintf(1,'Distance from object to the optical center:    so = %3.5f\n',so);
fprintf(1,'Distance from image to the optical center:     si = %3.5f\n',si);
fprintf(1,'Direction of the calibrated plane:      dir_plane = [%3.5f, %3.5f, %3.5f]\n',dir_plane);
fprintf(1,'Angle between the plane and the image:      theta = %3.5f in degree\n',theta);
fprintf(1,'--------------------------------------------------------------------------------\n');

