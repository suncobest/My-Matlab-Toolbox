% center = [-5 -5 -5 0 0 0 5 5 5 0;-4 0 4 -4 0 4 -4 0 4 0;0 0 0 0 0 0 0 0 0 3];
% voxels = subvoxcenter(center,2);

% center = rodrigues([1;2;3])* center;
% voxels  = subvoxcentercenter,3,[5;4;3],[1;2;3]);
% % [Vcube,face] = gen_cuboids(center,[5;4;3],[1;2;3]);
% 
% [sVcube,sface] = gen_cuboids(voxels,[5;4;3]/3,[1;2;3]);

n = 10;
center = 20 * randn(3,n);
nxyz = randi(4,3,n);
% nxyz = 1; snum = 1;
abc = randi(20,3,n);   % 1 row to make cubes; 3 row to extract cuboids; 
om = randn(3,n);
voxels = subvoxcenter(center,nxyz,abc,om);

snum = prod(nxyz,1);
sabc = abc./nxyz;
scube_abc = cell(1,n);
scube_om = scube_abc;
for kk=1:n,
    scube_abc{kk} = reshape(repmat(sabc(:,kk),1,snum(kk)),3,[]);
    scube_om{kk} = reshape(repmat(om(:,kk),1,snum(kk)),3,[]);
end;
scube_abc= cell2mat(scube_abc);
scube_om= cell2mat(scube_om);

[sVcube,sface] = gen_cuboids(voxels,scube_abc,scube_om);


% n = size(face,2);

n = size(sface,2);
cdata = jet(n);

figure(1), cameratoolbar, axis vis3d, hold on; grid on; 
plot3(center(1,:),center(2,:),center(3,:),'go','linewidth',2)
plot3(voxels(1,:),voxels(2,:),voxels(3,:),'b+')

% title('\bf\fontname{Times New Roman}\fontsize{16}\color{black}Draw cubes demo');


%%

% p = patch('Faces', face', 'Vertices', Vcube', 'FaceColor', 'b', 'EdgeColor','k','FaceAlpha',0.5);

p = patch('Faces', sface', 'Vertices', sVcube', 'FaceColor', 'b', 'EdgeColor','k','FaceAlpha',0.5);
light('Position',3*[0.5 -0.8 1],'Style','local'); 
material shiny 
camlight(20,0)  % camlight('headlight')
% light('Position',[1 -1 1],'Style','infinite');

set(p, 'FaceColor','flat','FaceVertexCData',cdata,'FaceAlpha',0.3)
% plot3(center(1,:),center(2,:),center(3,:),'r+')

% clear cdata 
% set(gca,'CLim',[0 50])  % 'CLim'º¥[cmin, cmax]£¨Color axis limits. £®Axis Ù–‘£©
% cdata = (1:n)';
% set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','direct')

axis equal tight, axis off;
