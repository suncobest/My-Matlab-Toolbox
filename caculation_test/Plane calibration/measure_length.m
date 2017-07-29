% Measurements of a line segment
fprintf(1,'\nMeasure a line segment:\nClick on any two points on the calibrated plane(vector direction: P1--->P2)...\n');
figure(2);
image(I);
colormap(gray(256));
set(2,'color',[1 1 1]);
title('Click on any two points on the calibrated plane(direction: P1--->P2)');
axis image;
hold on;
x = [];y = [];
for count = 1:2
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
pedge = [x(2)-x(1), y(2)-y(1),0];
lpedge = norm(pedge);

% dir_prj = dir_plane(1:2);
% dir_prj = [dir_prj/norm(dir_prj),0];
% pedge_dir = sum(pedge.*dir_prj);
% edge = pedge_dir; 
% pedge_dir = pedge_dir * dir_prj; 
% pedge_ndir = pedge - pedge_dir;
% ndir_prj = pedge_ndir/norm(pedge_ndir);
% hand = cross(dir_prj, ndir_prj);
% hand = hand(3);
% edge = (edge/dir_plane(3)* hand * cross(ndir_prj,dir_plane) + pedge_ndir)*pix/mag;

edge3d = invprj_vector(pedge,dir_plane) * pix /mag;
ledge3d = norm(edge3d);

if exist('Homo','var')
   Xw = Homo\([x,y,ones(2,1)]');
   Xw = Xw./(ones(3,1)*Xw(3,:));
   edge2d = Xw(1:2,2)-Xw(1:2,1);
   ledge2d = norm(edge2d);
end

fprintf(1,'\nMeasurement results:\n---------------------------------------------------------------------------------------\n');
fprintf(1,'Image length of the line segment in pixel:        lpedge = %3.5f\n',lpedge);
fprintf(1,'Physical length of the line segment(3D):         ledge3d = %3.5f\n',ledge3d);
fprintf(1,'3D vector of the line segment in the scene:       edge3d = [%3.5f, %3.5f, %3.5f]\n\n',edge3d);
if exist('Homo','var')
    fprintf(1,'From homography we have:\n');
    fprintf(1,'Physical length of the line segment(2D):         ledge2d = %3.5f\n',ledge2d);
    fprintf(1,'2D vector of the line segment on the plane:       edge2d = [%3.5f, %3.5f]\n',edge2d);
end
fprintf(1,'---------------------------------------------------------------------------------------\n');
figure(2);
if exist('oxy','var')
    delta = min(nx,ny)/30;
    arrow3(oxy(:,1)',oxy(:,2)','b2',1.2,2);
    arrow3(oxy(:,1)',oxy(:,4)','b2',1.2,2);
    vX = oxy(:,2)-oxy(:,1);
    vX = vX/norm(vX);
    vY = oxy(:,4)-oxy(:,1);
    vY = vY/norm(vY);
    vO = oxy(:,1)-oxy(:,3);
    vO = vO/norm(vO);
    position_x = oxy(:,2) + delta*vX;
    position_y = oxy(:,4) + delta*vY;
    position_o = oxy(:,1) + delta*vO;
    text('Interpreter','latex','String','$$x$$','Position',position_x','FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
    text('Interpreter','latex','String','$$y$$','Position',position_y','FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
    text('Interpreter','latex','String','$$o$$','Position',position_o','FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
end
if exist('ledge2d','var')
    text(sum(x)/2, sum(y)/2, num2str(ledge2d), 'color','g','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
else
    text(sum(x)/2, sum(y)/2, num2str(ledge3d), 'color','g','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
end
hold off;

