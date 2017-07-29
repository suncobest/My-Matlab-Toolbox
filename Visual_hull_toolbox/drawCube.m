center = [-5 -5 -5 0 0 0 5 5 5 0;-4 0 4 -4 0 4 -4 0 4 0;0 0 0 0 0 0 0 0 0 3];
% [Vcube,face] = gen_cuboids(center,[5;4;3]);

[Vcube,face] = gen_cuboids(center,3);

% center = rodrigues([1;2;3])* center;
% [Vcube,face] = gen_cuboids(center,[5;4;3],[1;2;3]);

% n = 10;
% center = 10 * randn(3,n);
% abc = randi(10,3,n);   % 1 row to make cubes; 3 row to extract cuboids; 
% om = randn(3,n);
% [Vcube,face] = gen_cuboids(center,abc,om);

n = size(face,2);
cdata = jet(n);

figure(1), cameratoolbar, axis vis3d, hold on; grid on; 
% title('\bf\fontname{Times New Roman}\fontsize{16}\color{black}Draw cubes demo');

p = patch('Faces', face', 'Vertices', Vcube', 'FaceColor', 'b', 'EdgeColor','k','FaceAlpha',0.8);
light('Position',3*[0.5 -0.8 1],'Style','local'); 
material shiny 
camlight(20,0)  % camlight('headlight')
% light('Position',[1 -1 1],'Style','infinite');

set(p, 'FaceColor','flat','FaceVertexCData',cdata,'FaceAlpha',0.7)
% plot3(center(1,:),center(2,:),center(3,:),'r+')

clear cdata 
set(gca,'CLim',[0 50])  % 'CLim'即[cmin, cmax]，Color axis limits. （Axis属性）
cdata = (1:n)';
set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','direct')


axis equal tight, axis off;

filename = 'cube demo3.gif';
view(3)
% [az,el]=view;
az0 = 0;
el = 20;
dt = 0.1;  % gif 的帧时间

%% --------------------------------------------------------------------------------------%
%%% 沿az方向往返旋转

daz = 2; % 0.8   % 每一帧方位角（azimuth）的增量
n = 20; % 3  % 单趟帧数

for i = 1 : n
    view(az0 + daz*(i-1),el)
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if  i == 1
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
    end
end

for i = n-1:-1:2
    view(az0 + daz*(i-1),el)
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
	imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
end

%% --------------------------------------------------------------------------------------%
%%% 沿az方向旋转一周

% daz = 5;  % 每一帧方位角（azimuth）的增量
% n = 360/daz;   % 旋转一周步数
% 
% for i = 1 : n-1
%     view(az0 + daz*(i-1),el)
%     drawnow
%     frame = getframe(1);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
%     if  i == 1
% 		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',dt);
% 	else
% 		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',dt);
%     end
% end

% print(1,'-dpng','cube demo.png');

