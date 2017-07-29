% Measurements of two connected line segments and the angle between them
fprintf(1,'\nMeasure two connected line segments and the angle between them:\n');
figure(2);
image(I);
colormap(gray(256));
set(2,'color',[1 1 1]);
title('Click on the three vertexes of an angle in the counter-clockwise direction');
disp('Click on the three vertexes of an angle in the counter-clockwise direction on the calibrated plane(the second vertex is the corner)...');
axis image;
hold on;
x= [];y = [];
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

pedge1 = [x(1)-x(2), y(1)-y(2),0];
pedge2 = [x(3)-x(2), y(3)-y(2),0];
lpedge1 = norm(pedge1);
lpedge2 = norm(pedge2);
pangle1_2 = acosd(sum(pedge1.*pedge2)/(lpedge1*lpedge2));
edge1 = invprj_vector(pedge1,dir_plane) * pix/mag;
edge2 = invprj_vector(pedge2,dir_plane) * pix/mag;
ledge1 = norm(edge1);
ledge2 = norm(edge2);
angle1_2 = acosd(sum(edge1.*edge2)/(ledge1*ledge2));
    
if exist('Homo','var')
    Xw = Homo\([x,y,ones(3,1)]');
    Xw = Xw./(ones(3,1)*Xw(3,:));
    l1 = Xw(1:2,1)-Xw(1:2,2);
    l2 = Xw(1:2,3)-Xw(1:2,2);
    ll1 = norm(l1);
    ll2 = norm(l2);
    angle12 = acosd(sum(l1.*l2)/(ll1*ll2));
end

fprintf(1,'\nMeasurement results:\n---------------------------------------------------------------------------------------\n');
disp('Parameters of the image angle:');
fprintf(1,'Image length of the edge1 in pixel:             lpedge1 = %3.5f\n',lpedge1);
fprintf(1,'Image length of the edge2 in pixel:             lpedge2 = %3.5f\n',lpedge2);
fprintf(1,'Image angle between edge1 and edge2:          pangle1_2 = %3.5f in degree\n\n',pangle1_2);
disp('Parameters of the physical angle');
fprintf(1,'Physical length of the edge1(3D):                ledge1 = %3.5f\n',ledge1);
fprintf(1,'Physical length of the edge2(3D):                ledge2 = %3.5f\n',ledge2);
fprintf(1,'Physical angle between edge1 and edge2(3D):    angle1_2 = %3.5f in degree\n\n',angle1_2);
if exist('Homo','var')
    fprintf(1,'From homography we have:\n');
    fprintf(1,'Physical length of the edge1(2D):                   ll1 = %3.5f\n',ll1);
    fprintf(1,'Physical length of the edge2(2D):                   ll2 = %3.5f\n',ll2);
    fprintf(1,'Physical angle between edge1 and edge2(2D):     angle12 = %3.5f\n',angle12);
end
fprintf(1,'---------------------------------------------------------------------------------------\n');

delta = min(nx,ny)/30;
lu_corner = delta*[1,1];
radius = 3*delta;

dstart = atan2d(pedge1(2),pedge1(1));
dend = atan2d(pedge2(2),pedge2(1));

if dend-dstart>180,
    dend=dend-360;
end;

position_angle = [x(2),y(2)] + (delta+radius)*[cosd((dstart+dend)/2), sind((dstart+dend)/2)];

AX = figure(2);
image(I);
colormap(gray(256));
set(AX,'color',[1 1 1]);
title('Measurement results of the physical angle');
axis image;
hold on;

if 0 && exist('oxy','var')
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

draw_arc(AX,[x(2),y(2)],radius,dstart,dend,20);
arrow3([x(2),y(2)],[x(1),y(1)],'b2',1.2,2);
arrow3([x(2),y(2)],[x(3),y(3)],'b2',1.2,2);
position_l1 = [x(1),y(1)] + delta*pedge1(1:2)/lpedge1;
position_l2 = [x(3),y(3)] + delta*pedge2(1:2)/lpedge2;

text('Interpreter','latex','String','$$l_1$$','Position',position_l1,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
text('Interpreter','latex','String','$$l_2$$','Position',position_l2,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');
text('Interpreter','latex','String','$$\alpha$$','Position',position_angle,'FontSize',20,'color',[0 1 0],'fontweight','bold','HorizontalAlignment','center');

if exist('Homo','var')
    text('Interpreter','latex','String',['$$l_1 = ' num2str(ll1) '$$'],'Position',lu_corner, 'color','g','Fontsize',20,'fontweight','bold');
    text('Interpreter','latex','String',['$$l_2 = ' num2str(ll2) '$$'],'Position',lu_corner+[0 2*delta], 'color','g','Fontsize',20,'fontweight','bold');
    text('Interpreter','latex','String',['$$\alpha = ' num2str(angle12) '$$'],'Position',lu_corner+[0 4*delta], 'color','g','Fontsize',20,'fontweight','bold');
    
else
    text('Interpreter','latex','String',['$$l_1 = ' num2str(ledge1) '$$'],'Position',lu_corner, 'color','g','Fontsize',20,'fontweight','bold');
    text('Interpreter','latex','String',['$$l_2 = ' num2str(ledge2) '$$'],'Position',lu_corner+[0 2*delta], 'color','g','Fontsize',20,'fontweight','bold');
    text('Interpreter','latex','String',['$$\alpha = ' num2str(angle1_2) '$$'],'Position',lu_corner+[0 4*delta], 'color','g','Fontsize',20,'fontweight','bold');
end

hold off;