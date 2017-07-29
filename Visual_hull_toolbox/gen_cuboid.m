function [vertices,polyface] = gen_cuboid(center,abcsize,omlist)

% center,abcsize,omlist ----- (3,N)
% vertices ------- (3,8*N)
% polyface ------- (4,6*N)
% N is the total number of cuboids. center are the center of every cuboids. 
% abcsize are the sizes of cuboids. omlist are the directions of cuboids.
% Every cuboid has 8 vertices and 6 faces(polygon)

cube = [0 1 1 0 0 1 1 0
        0 0 1 1 0 0 1 1
        0 0 0 0 1 1 1 1]-ones(3,8)/2;    % unit cube vertices initialization(1-2-3-4_5-6-7-8)
    
face = [1  5  1  4  2  1
        4  6  2  8  3  5
        3  7  6  7  7  8
        2  8  5  3  6  4];   %  six faces with out-pointed normal vector （right-hand rule）
    
% face的元素为顶点的指标索引，每一列表示一个多边形，按顺序连接顶点构成边，顶点顺序可以滚动或反向，但不可交叉    
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

polyface = zeros(4,6*np);   % 所有长方体的6个面（四个顶点的指标）
for i = 1:np,
    polyface(:,6*(i-1)+1:6*i) = 8*(i-1)+face;  % generate faces (idx of vertices)  
end;

if nargin==1,
    vertices = repmat(cube,1,np) + Tlist;  % translation
    return;
end;

[m,n] = size(abcsize);
if m==1,
    abcsize = abcsize(ones(3,1),:);
elseif m ~= 3,
     error('Size of cubes have to be 1 or 3 dimensional!');
end;

% resize cubes
if n==1;
    vertices = diag(abcsize) * repmat(cube,1,np);  % all cubes have the same size
else
    if n~=np,
        error('Points number unmatched!');
    end;
    for i = 1:np,
        vertices(:,8*(i-1)+1:8*i) = diag(abcsize(:,i)) * cube;  % every cube has its own size
    end;
end;

if nargin == 2,
    vertices = vertices + Tlist;  % translation
    return;
end;

% rotate cubes
[m,n] = size(omlist);
if m ~= 3,
    error('Omega (rotation) have to be 3 dimensional! (3 rows)');
end;

if n ==1,
    vertices = rodrigues(omlist) * vertices;
else
    if n~=np,
        error('Points number unmatched!');
    end;
    for i = 1:np,
        vertices(:,8*(i-1)+1:8*i) =  rodrigues(omlist(:,i)) * vertices(:,8*(i-1)+1:8*i);
    end;
end;

% translate cubes
vertices = vertices + Tlist;

return;

