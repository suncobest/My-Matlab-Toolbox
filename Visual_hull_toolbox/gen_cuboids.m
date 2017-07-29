function [vertices,polyface] = gen_cuboids(center,abcsize,omlist)

% center,abcsize,omlist ----- (3,N)
% vertices ------- (3,8*N)
% polyface ------- (4,6*N)
% N is the total number of cuboids. center are the center of every cuboids.
% abcsize are the sizes of cuboids. omlist are the directions of cuboids.
% Every cuboid has 8 vertices and 6 faces(polygon)
% vertices are the 8 vertex of every cuboid, polyfaces are the 6 faces of
% every cuboid. Each face contain column index of its four vertices.

cube = [0, 1, 1, 0, 0, 1, 1, 0;
        0, 0, 1, 1, 0, 0, 1, 1;
        0, 0, 0, 0, 1, 1, 1, 1]-ones(3,8)/2;    % unit cube vertices initialization(1-2-3-4_5-6-7-8)

face = [1, 5, 1, 4, 2, 1;
        4, 6, 2, 8, 3, 5;
        3, 7, 6, 7, 7, 8;
        2, 8, 5, 3, 6, 4];   %  six faces with out-pointed normal vector (right-hand rule)

% face的元素为顶点的指标索引，每一列表示一个多边形，
% 按顺序连接顶点构成边，顶点顺序可以滚动或反向，但不可交叉    
% patch('Faces', face', 'Vertices', cube', 'FaceColor', 'b', 'EdgeColor', 'k');  

if nargin==0,
    vertices = cube;
    polyface = face;
    return;
end;

[m,np] = size(center);
if m ~= 3,
    error('Points have to be 3 dimensional! (3 rows)');
end;

Tlist = reshape(repmat(center,8,1),3,[]);  % 平移8个顶点需要将T复制8遍

if nargout==2,
    polyface = reshape(repmat(face(:),1,np)+repmat(8*(0:np-1),24,1), 4, []); % generate faces (idx of vertices)
end;

if nargin==1,
    vertices = repmat(cube,1,np) + Tlist;  % translation
    return;
end;

% resize cubes  (make vertices an array size of [3,8,N])
[m,n] = size(abcsize);

if m==1,
    if n==1,
        verts = abcsize * repmat(cube,[1,1,np]);  % all cubes have the same size 
    elseif n==np,
        verts = repmat(cube,[1,1,np]) .* repmat(reshape(abcsize,1,1,[]),3,8);   % every cube has its own size
    else
        error('Points number unmatched!');
    end;
elseif m==3,
    if n==1,
        verts = repmat((repmat(abcsize,1,8) .* cube),[1,1,np]);  % all cubes have the same size
    elseif n==np,
        verts = repmat(cube,[1,1,np]) .* repmat(reshape(abcsize,3,1,[]),1,8);  % every cube has its own size 
    else
        error('Points number unmatched!');
    end;
else
    error('Size of cubes have to be 1 or 3 dimensional!');
end;


if nargin == 2,
    vertices = reshape(verts,3,[]) + Tlist;  % translation
    return;
end;

% rotate cubes
[m,n] = size(omlist);
if m ~= 3,
    error('Omega (rotation) have to be 3 dimensional! (3 rows)');
end;

% all cubes have the same orientation (对于计算一个om，axis2mat也可以换成rodrigues)
if n ==1,
    vertices = rodrigues(omlist) * reshape(verts,3,[]);   
else
    if n==np,
        vertices = zeros(3,8*np);
        for i = 1:np,
            vertices(:,8*(i-1)+1:8*i) =  rodrigues(omlist(:,i)) * verts(:,:,i);
        end;
    else
        error('Points number unmatched!');
    end;
end;

% translate cubes
vertices = vertices + Tlist;

return;


%% Test
n = 3;
center = 10 * randn(3,n);
abc = randi(10,3,n);   % 1 row to make cubes; 3 row to extract cuboids; 
om = randn(3,n);
tic;
[Vcubes,faces] = gen_cuboids(center,abc,om);
t1 = toc
% clear Vcubes faces
tic;
[Vcube,face] = gen_cuboid(center,abc,om);
t2 = toc
